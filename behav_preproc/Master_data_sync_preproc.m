function Master_data_sync_preproc(sessions, do_plot, project_hint)
% Master_data_sync_preproc
% Unified wrapper for React_Passive (RP), React_Active (RA) and Tonotopy (T): sync all data streams and build epochs.
%   - sync triggers (csv and/or OE)
%   - parse baphy m-file
%   - regularize TTLs
%   - build epochs
%
% INPUT
%   sessions : cell array of session folders, e.g.
%       sessions = {
%         '/home/Arsenii/React_Passive/Processed_data/Edel/20220419_1_m_C'
%         '/home/Arsenii/React_Passive/Processed_data/Kosichka/20220309_1_m_A'
%         '/home/Arsenii/React_Passive/Processed_data/Chabichou/20220301_1_m_A'
%         ...
%       };
close all
if nargin < 2 || isempty(do_plot)
    do_plot = false;
end

if nargin < 3
    project_hint = '';
end

if ~iscell(sessions)
    error('Input "sessions" must be a cell array of session paths');
end

for sess = 1:numel(sessions)
    datapath = sessions{sess};
    disp('------------------------------------------')
    disp(['Working on ' datapath])
    
    if ~isfolder(datapath)
        warning(['Folder not found, skipping: ' datapath])
        continue
    end
    
    try
        % 1) detect project
        project = detect_project_from_path(datapath, project_hint);
        fprintf(' [1] Project = %s\n', project);
        
        % 2) Check which modalities are available
        lfp_dir = fullfile(datapath,'ephys','LFPData');
        if ~isfolder(lfp_dir)
            lfp_dir = fullfile(datapath,'LFPData'); % fallback for older layout
        end
        has_OE = isfolder(lfp_dir) && ~isempty(dir(fullfile(lfp_dir,'LFP*.mat')));
        
        % 3) Passive-specific csv sync for Edel and Chabichou (and other csv-only cases which probably will never exist)
        trig_csv = [];
        if strcmp(project,'RP')
            if contains(datapath,'Edel') || contains(datapath,'Chabichou')
                disp('  [3] Syncing Baphy�fUS from csv (RP_sync_triggers_passive)...')
                trig_csv = RP_sync_triggers_passive(datapath);
            else
                disp('  [3] No csv sync step for this animal (skip RP_sync_triggers_passive).')
            end
        end
        
        % 4-5) Get channel configuration and extract TTLs from OE for this session
        trigOE = [];
        if has_OE
            disp('  [4] Loading trigger channel config (get_trigger_config)...')
            cfg = get_trigger_config(datapath);
            cfg.plt = do_plot;
            cfg.project_hint = project_hint;
            
            disp('  [5] Extracting OE TTL triggers (extract_triggers_oe)...')
            trigOE = extract_triggers_oe(datapath, cfg);
        else
            disp('  [4-5] No LFPData/OpenEphys for this session, skipping OE TTL extraction.');
            cfg.project_hint = project_hint;
        end
        
        % 6) Parse Baphy m-file for React Passive (categories, timings, etc.)
        if strcmp(project,'RP')
            disp('  [6] Parsing Baphy file (RP_parse_baphy_passive)...')
            Baphy = RP_parse_baphy_passive(datapath);
        elseif strcmp(project,'RA')
            disp('  [6] Parsing Baphy file (RA_parse_baphy_active)...')
            Baphy = RA_parse_baphy_active(datapath, trigOE);
        elseif strcmp(project,'T')
            disp('  [6] Parsing Baphy file (RP_parse_baphy_tonotopy)...')
            Baphy = RP_parse_baphy_tonotopy(datapath);
        end
    
        % 6b) Regularise / check TTLs
        if has_OE && ~isempty(trigOE)
            disp('  [6b] Regularising OE TTLs (regularize_all_ttl)...')
            try
                trigOE = regularize_all_ttl(trigOE, Baphy, datapath);
            catch MEreg
                warning('regularize_all_ttl crashed: %s', MEreg.message);
            end
        elseif ~has_OE && ~isempty(trig_csv)
            disp('  [6b] Regularising csv TTLs (regularize_all_ttl_csv)...')
            trig_csv = regularize_all_ttl_csv(trig_csv, Baphy, datapath);
        else
            disp('  [6b] No TTLs available for regularisation/checking.');
        end
        
        % 6c) Edel-specific: derive Baphy trial info from OE TTL bursts
        if has_OE && contains(datapath,'Edel')
            disp('  [6c] Extracting Edel Baphy trial info from OE TTLs...')
            trigOE = RP_add_baphy_trial_info_from_ttl(trigOE, trig_csv, Baphy, datapath);
        end
        
        % 7) Build epochs on the OE/common time axis
        Epochs = [];
        if strcmp(project,'RP')
            if has_OE && ~isempty(trigOE)
                disp('  [7] Building epochs (RP_build_epochs_passive)...')
                Epochs = RP_build_epochs_passive(datapath, trigOE, Baphy);
            else
                disp('  [7] Building epochs from csv only (RP_build_epochs_passive_csv_only)...')
                if isempty(trig_csv)
                    error('Master_data_sync_preproc:NoCSVForCsvOnly', ...
                        'RP csv-only session but trig_csv is empty.');
                end
                Epochs = RP_build_epochs_passive_csv_only(datapath, trig_csv, Baphy);
            end
        elseif srcmp(project, 'Tonotopy')
            if has_OE && ~isempty(trigOE)
                disp('  [7] Building epochs (RP_build_epochs_tonotopy)...')
                Epochs = RP_build_epochs_tonotopy(datapath, trigOE, Baphy);
            else
                disp(' No tonotopy yet for .csv triggers')
            end
            
        else
            if ~has_OE || isempty(trigOE)
                error('RA requires OE triggers. No OE LFPData detected.');
            end
            disp('  [7] Building epochs (RA_build_epochs_active)...')
            Epochs = RA_build_epochs_active(datapath, trigOE, Baphy);
        end
        
        % 8) Optional trial-count sanity check (now using scalar counts)
        if exist('check_trial_counts','file') == 2
            try
                disp('[8] Checking trial counts (check_trial_counts)...')
                
                % n_oe: from OE Baphy TTLs, if available
                n_oe = [];
                
                if ~isempty(trigOE) && isfield(trigOE,'baphy')
                    if isfield(trigOE.baphy,'t_s') && ~isempty(trigOE.baphy.t_s)
                        n_oe = numel(trigOE.baphy.t_s);
                    elseif isfield(trigOE.baphy,'t_raw_s') && ~isempty(trigOE.baphy.t_raw_s)
                        [t_on_oe, ~] = detect_trial_onsets_from_baphy_ttl(datapath, trigOE.baphy.t_raw_s, Baphy.n_trials);
                        n_oe = numel(t_on_oe);
                    elseif isfield(trigOE.baphy,'trial_start_ts') && ~isempty(trigOE.baphy.trial_start_ts)
                        n_oe = numel(trigOE.baphy.trial_start_ts);
                    end
                end
                
                % n_csv: from csv Baphy TTLs, if available
                n_csv = [];

                if ~isempty(trig_csv) && isfield(trig_csv,'baphy_time_s') && ~isempty(trig_csv.baphy_time_s)
                    [t_on_csv, ~] = detect_trial_onsets_from_baphy_ttl(datapath, trig_csv.baphy_time_s, Baphy.n_trials);
                    n_csv = numel(t_on_csv);
                end
                
                % n_baphy: from Baphy metadata
                n_baphy = Baphy.n_trials;
                
                check_trial_counts(n_oe, n_baphy, n_csv, datapath, 'RP_passive');
                
            catch ME_chk
                warning('check_trial_counts crashed: %s', ME_chk.message);
            end
        else
            disp('  [4e] check_trial_counts.m not found, skipping.');
        end
        
        % Small textual summary
        n_csv_baphy = 0;
        if ~isempty(trig_csv) && isfield(trig_csv,'baphy_time_s')
            n_csv_baphy = numel(trig_csv.baphy_time_s);
        end
        n_oe_baphy = 0;
        if ~isempty(trigOE) && isfield(trigOE,'baphy') && isfield(trigOE.baphy,'t_raw_s')
            n_oe_baphy = numel(trigOE.baphy.t_raw_s);
        end
        fprintf('  Summary: Baphy-trials=%d | csv-Baphy-TTL=%d | OE-Baphy-TTL=%d\n', ...
            Baphy.n_trials, n_csv_baphy, n_oe_baphy);
        
        
        % 6) Save everything in one place
        out_file = fullfile(datapath, 'Master_sync.mat');
        save(out_file, 'trig_csv', 'trigOE', 'Baphy', 'Epochs');
        disp(['  Saved master sync file: ' out_file])
        
    catch ME
        warning('Error in session %s:\n  %s', datapath, ME.message);
        disp(getReport(ME,'basic'))
        continue
    end
end

disp('Master_data_sync_preproc finished.')

end

%             if do_plot && isfield(trigOE,'fus') && isfield(trigOE.fus,'t_s') ...
%                     && ~isempty(trigOE.fus.t_s)
%                 figure('Name',['OE fUS raw vs regularised TTLs - ' datapath],'NumberTitle','off');
%                 hold on 
%                 plot(trigOE.fus.t_s, 1.1, 'r*');
%                 plot(trigOE.fus.t_raw_ts/1e4, 1, 'b*');
%                 plot(trigOE.fus.t_s(find(ismember(trigOE.fus.t_s, trigOE.fus.t_raw_ts/1e4)==0)), 0.99, 'g*');
%                 xlabel('Time (s)'); ylim([0 2]);
%                 title('OE fUS raw vs regularised TTLs');
%             end
%             
%                         
%             if do_plot && isfield(trig_csv,'fus_time_s') && ~isempty(trig_csv.fus_time_s)
%                 figure('Name',['csv fUS TTL - ' datapath],'NumberTitle','off');                
%                 plot(trig_csv.fus_time_s, ones(size(trig_csv.fus_time_s)), 'k*');
%                 xlabel('Time (s)'); ylim([0.5 1.5]);
%                 title('csv fUS TTLs (after regularisation)');
%             end