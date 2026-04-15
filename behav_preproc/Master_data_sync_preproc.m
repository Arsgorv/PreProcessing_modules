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
                disp('  [3] Syncing Baphy and fUS from csv (RP_sync_triggers_passive)...')
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
            % RA training (single mfile) vs RA experiment (Conditioning+PostTest)
            isRAexp = ~isempty(dir(fullfile(datapath,'ephys','*_RA_Conditioning*'))) || ...
                ~isempty(dir(fullfile(datapath,'ephys','*_RA_PostTest*')));
            
            if ~isRAexp
                stimdir = fullfile(datapath,'stim');
                mAll = dir(fullfile(stimdir,'*.m'));
                if numel(mAll) <= 1
                    disp('  [6] Parsing Baphy file (RA_parse_baphy_active)...')
                    Baphy = RA_parse_baphy_active(datapath, trigOE);
                    disp('  [7] Building epochs (RA_build_epochs_active)...')
                    Epochs = RA_build_epochs_active(datapath, trigOE, Baphy);
                else
                    disp('  [6] RA training: multiple mfiles detected, parsing/stitching conservatively...')
                    [Baphy, Epochs, partsTrain] = RA_train_parse_multimfile(datapath, trigOE, 20, 5);
                    save(fullfile(datapath,'Master_sync_training_parts.mat'), 'partsTrain');
                end
            else
                disp('  [6] RA experiment: building run manifest + parsing Conditioning/PostTest separately')
                
                mf = fullfile(datapath,'analysis','run_manifest_RAexp.mat');
                if exist(mf,'file')
                    Srun = load(mf,'RUN'); RUN = Srun.RUN;
                else
                    RUN = RAE_make_run_manifest(datapath);
                end
                
                Epochs = struct();
                Epochs.highlevel = struct();
                phList = {'PreSleep','Conditioning','PostSleep','PostTest'};
                for p = 1:numel(phList)
                    ph = phList{p};
                    t0 = RUN.run.(ph).t0_ts; t1 = RUN.run.(ph).t1_ts;
                    if isfinite(t0) && isfinite(t1)
                        Epochs.highlevel.(ph) = intervalSet(t0, t1);
                    end
                end
                
                % ---- Conditioning ----
                stimdir = fullfile(datapath,'stim');
                mCond = dir(fullfile(stimdir,'*_Conditioning*.m'));
                mPost = dir(fullfile(stimdir,'*_PostTest*.m'));
                if ~isempty(mCond)
                    [BC, EC, trigC, partsC] = RAExp_parse_phase(datapath, trigOE, RUN, 'Conditioning', ...
                        'Conditioning', 20, 5);  % blockGap_s=20, pad_s=5
                    
                    if ~isempty(BC)
                        BC = RA_fix_baphy_trial_masks(BC);
                        Epochs.Conditioning = EC;
                        save(fullfile(datapath,'Master_sync_Conditioning.mat'), 'trigC','BC','EC','partsC');
                    end
                end
                
                % ---- PostTest ----
                if ~isempty(mPost)
                    [BP, EP, trigP, partsP] = RAExp_parse_phase(datapath, trigOE, RUN, 'PostTest', ...
                        'PostTest', 20, 5);
                    
                    if ~isempty(BP)
                        BP = RA_fix_baphy_trial_masks(BP);
                        Epochs.PostTest = EP;
                        save(fullfile(datapath,'Master_sync_PostTest.mat'), 'trigP','BP','EP','partsP');
                    end
                end
                
                % Optional: build union epochs for convenience
                if isfield(Epochs,'Conditioning') && isfield(Epochs,'PostTest')
                    Epochs.trial_all = union(Epochs.Conditioning.trial_all, Epochs.PostTest.trial_all);
                    Epochs.stim_all  = union(Epochs.Conditioning.stim_all,  Epochs.PostTest.stim_all);
                end
                
                % store a compact Baphy container (keeps compatibility: you still have per-phase mats)
                Baphy = struct();
                if exist('BC','var')
                    BC = RA_fix_baphy_trial_masks(BC);
                    Baphy.Conditioning = BC; 
                end
                if exist('BP','var')
                    BP = RA_fix_baphy_trial_masks(BP);
                    Baphy.PostTest     = BP;
                end
                
                % Provide a legacy-compatible total trial count (used by regularize/checks/summary)
                Baphy.n_trials = 0;
                if isfield(Baphy,'Conditioning') && isfield(Baphy.Conditioning,'n_trials')
                    Baphy.n_trials = Baphy.n_trials + Baphy.Conditioning.n_trials;
                end
                if isfield(Baphy,'PostTest') && isfield(Baphy.PostTest,'n_trials')
                    Baphy.n_trials = Baphy.n_trials + Baphy.PostTest.n_trials;
                end
                
                Baphy.isRAexp = true;
            end
            
            
            
        elseif strcmp(project,'T')
            disp('  [6] Parsing Baphy file (RP_parse_baphy_tonotopy)...')
            Baphy = RP_parse_baphy_tonotopy(datapath);
        end
        
        % 6b) Regularise / check TTLs
        if has_OE && ~isempty(trigOE)
            disp('  [6b] Regularising OE TTLs (regularize_all_ttl)...')
            if strcmp(project,'Arousal')
                trigOE = regularize_all_ttl_Arousal(trigOE, datapath);
            else
                try
                    trigOE = regularize_all_ttl(trigOE, Baphy, datapath);
                catch MEreg
                    warning('regularize_all_ttl crashed: %s', MEreg.message);
                end
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
        elseif strcmp(project, 'Tonotopy')
            if has_OE && ~isempty(trigOE)
                disp('  [7] Building epochs (RP_build_epochs_tonotopy)...')
                Epochs = RP_build_epochs_tonotopy(datapath, trigOE, Baphy);
            else
                disp(' No tonotopy yet for .csv triggers')
            end
        elseif strcmp(project, 'Arousal')
            if has_OE && ~isempty(trigOE)
                disp('  [7] Building epochs (RP_build_epochs_Arousal)...')
                Epochs = RP_build_epochs_arousal(datapath, trigOE);
            else
                disp(' No Arousal sync yet for .csv triggers')
            end
            
        else
            % RA: epochs already built in the RA branch above (training OR experiment)
            disp('  [7] RA: epochs already built, skipping.');
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
                n_baphy = get_baphy_n_trials_any(Baphy);
                
                check_trial_counts(n_oe, n_baphy, n_csv, datapath, 'RP_passive');
                
            catch ME_chk
                warning('check_trial_counts crashed: %s', ME_chk.message);
            end
        else
            disp('  [4e] check_trial_counts.m not found, skipping.');
        end
        
        if ~strcmp(project, 'Arousal')
            % Small textual summary baphy
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
        end 
        
        
        % 9) Save everything in one place
        out_file = fullfile(datapath, 'Master_sync.mat');
        if strcmp(project, 'Arousal')
            save(out_file, 'trig_csv', 'trigOE', 'Epochs');
        else
            save(out_file, 'trig_csv', 'trigOE', 'Baphy', 'Epochs');
        end
        disp(['  Saved master sync file: ' out_file])
        
    catch ME
        warning('Error in session %s:\n  %s', datapath, ME.message);
        disp(getReport(ME,'basic'))
        continue
    end
    close all
end

disp('Master_data_sync_preproc finished.')

end
% helper: subset trigOE fields by run time window
% (ticks are 1e4-compatible with  LFP time base)
function tr = subset_trigOE_to_window(trigOE, t0_ts, t1_ts)
% Subset trigOE modalities to [t0_ts, t1_ts] in ticks (1e4).
% Works if modalities store times in t_raw_ts and/or t_raw_s and/or trial_start_ts.

TsRate = 1e4;
t0_s = t0_ts / TsRate;
t1_s = t1_ts / TsRate;

tr = trigOE;
fns = {'baphy','face_cam','eye_cam','respi','heart','OneBox','fus'};

for ii = 1:numel(fns)
    fn = fns{ii};
    if ~isfield(tr,fn) || ~isstruct(tr.(fn)), continue; end
    
    S = tr.(fn);
    
    % --- subset t_raw_ts / t_raw_s consistently ---
    if isfield(S,'t_raw_ts') && ~isempty(S.t_raw_ts)
        x = S.t_raw_ts(:);
        keep = (x>=t0_ts & x<=t1_ts);
        S.t_raw_ts = x(keep);
        if isfield(S,'t_raw_s') && ~isempty(S.t_raw_s)
            xs = S.t_raw_s(:);
            if numel(xs) == numel(x)
                S.t_raw_s = xs(keep);
            else
                % fallback recompute from ts
                S.t_raw_s = S.t_raw_ts / TsRate;
            end
        else
            S.t_raw_s = S.t_raw_ts / TsRate;
        end
        
    elseif isfield(S,'t_raw_s') && ~isempty(S.t_raw_s)
        xs = S.t_raw_s(:);
        keep = (xs>=t0_s & xs<=t1_s);
        S.t_raw_s = xs(keep);
        if isfield(S,'t_raw_ts') && ~isempty(S.t_raw_ts)
            % if inconsistent, overwrite
            S.t_raw_ts = round(S.t_raw_s * TsRate);
        else
            S.t_raw_ts = round(S.t_raw_s * TsRate);
        end
    end
    
    % --- subset regularized t_s if present ---
    if isfield(S,'t_s') && ~isempty(S.t_s)
        xs = S.t_s(:);
        keep = (xs>=t0_s & xs<=t1_s);
        S.t_s = xs(keep);
    end
    
    % --- subset trial_start_ts / trial_stop_ts if present ---
    if isfield(S,'trial_start_ts') && ~isempty(S.trial_start_ts)
        y = S.trial_start_ts(:);
        kk = (y>=t0_ts & y<=t1_ts);
        S.trial_start_ts = y(kk);
    end
    if isfield(S,'trial_stop_ts') && ~isempty(S.trial_stop_ts)
        y = S.trial_stop_ts(:);
        kk = (y>=t0_ts & y<=t1_ts);
        S.trial_stop_ts = y(kk);
    end
    
    % --- update n_trials / n_pulses if applicable ---
    if isfield(S,'t_raw_ts') && ~isempty(S.t_raw_ts)
        if isfield(S,'n_trials')
            S.n_trials = numel(S.t_raw_ts);
        end
        if isfield(S,'n_pulses')
            S.n_pulses = numel(S.t_raw_ts);
        end
        if isfield(S,'n_raw_pulses')
            S.n_raw_pulses = numel(S.t_raw_ts);
        end
    end
    
    tr.(fn) = S;
end
end


function [Btrain, Etrain, parts] = RA_train_parse_multimfile(datapath, trigOE, blockGap_s, pad_s)
% Conservative training parser for sessions with multiple Baphy .m files.
% Updated v2:
%   - do not split training TTLs into gap-defined blocks
%   - consume the pruned trial TTL sequence in session order
%   - assign TTLs to each m-file sequentially using expected trial counts

Btrain = [];
Etrain = [];
parts = struct('mfile',{},'t0_s',{},'t1_s',{},'B',{},'E',{}, ...
    'n_trials_expected',{},'n_ttl_used',{});

stimdir = fullfile(datapath,'stim');
mfiles = RA_list_training_mfiles(stimdir);
if isempty(mfiles)
    warning('No training mfiles found in %s', stimdir);
    return
end

t = get_baphy_ttl_seconds(trigOE);
if isempty(t)
    warning('No Baphy TTLs found for multi-mfile training parse.');
    return
end

t = RA_prune_close_pulses(t, 0.5);
nTTL = numel(t);
nMf = numel(mfiles);

info = cell(nMf,1);
nExp = nan(nMf,1);
for k = 1:nMf
    info{k} = RA_get_mfile_info(mfiles{k});
    nExp(k) = info{k}.n_trials;
end

if all(isfinite(nExp))
    fprintf('  [RA training] multiple mfiles: expected total=%d, TTL=%d\n', sum(nExp), nTTL);
else
    fprintf('  [RA training] multiple mfiles: TTL=%d, some expected counts unavailable\n', nTTL);
end

Bparts = cell(nMf,1);
Eparts = cell(nMf,1);
cursor = 1;

for k = 1:nMf
    mf = mfiles{k};
    infoMf = info{k};

    if cursor > nTTL
        warning('[RA training] No TTLs left for mfile %d/%d: %s', k, nMf, mf);
        break
    end

    if isfinite(infoMf.n_trials) && infoMf.n_trials > 0
        stopIdx = min(nTTL, cursor + infoMf.n_trials - 1);
    else
        if k < nMf
            nLeftMf = nMf - k + 1;
            approxN = floor((nTTL - cursor + 1) / nLeftMf);
            stopIdx = min(nTTL, cursor + max(1, approxN) - 1);
        else
            stopIdx = nTTL;
        end
    end

    selT = t(cursor:stopIdx);

    if isempty(selT)
        warning('[RA training] Empty TTL assignment for mfile %d/%d: %s', k, nMf, mf);
        break
    end

    if isfinite(infoMf.n_trials) && numel(selT) < infoMf.n_trials
        warning('[RA training] mfile %d/%d expected %d trials but only %d TTLs remain.', ...
            k, nMf, infoMf.n_trials, numel(selT));
    end

    t0_s = selT(1) - pad_s;
    t1_s = selT(end) + pad_s;

    trK = subset_trigOE_to_window(trigOE, round(t0_s*1e4), round(t1_s*1e4));
    trK = RA_overwrite_baphy_ttls(trK, selT);

    outTag = sprintf('train_%d', k);
    Bk = RA_parse_baphy_active(datapath, trK, 'MFile', mf, 'OutTag', outTag);
    Ek = RA_build_epochs_active(datapath, trK, Bk);

    parts(k).mfile = mf;
    parts(k).t0_s = t0_s;
    parts(k).t1_s = t1_s;
    parts(k).B = Bk;
    parts(k).E = Ek;
    parts(k).n_trials_expected = infoMf.n_trials;
    parts(k).n_ttl_used = numel(selT);

    Bparts{k} = Bk;
    Eparts{k} = Ek;

    cursor = stopIdx + 1;
end

Bparts = Bparts(~cellfun(@isempty, Bparts));
Eparts = Eparts(~cellfun(@isempty, Eparts));
parts = parts(~cellfun(@isempty, {parts.mfile}));

if isempty(Bparts)
    warning('Could not parse any training multi-mfile parts.');
    return
end

if cursor <= nTTL
    warning('[RA training] %d trailing Baphy TTLs were left unused after sequential mfile assignment.', nTTL - cursor + 1);
end

Btrain = RAExp_stitch_baphy(Bparts);
Etrain = RAExp_stitch_epochs(Eparts);
end

function [Bphase, Ephase, trigPhase, parts] = RAExp_parse_phase(datapath, trigOE, RUN, phaseName, dlcSubdir, blockGap_s, pad_s)
% Parse a head-fixed phase that may have multiple Baphy mfiles and TTL blocks.
% Conservative updates:
%   - restrict to manifest phase window
%   - prune unrealistically close duplicate TTL pulses
%   - if a block has too many pulses, select the subset closest to the mfile

Bphase = [];
Ephase = [];
trigPhase = [];
parts = struct('mfile',{},'t0_s',{},'t1_s',{},'B',{},'E',{}, ...
    'n_trials_expected',{},'n_ttl_used',{});

TsRate = 1e4;
stimdir = fullfile(datapath,'stim');

mfiles = RAExp_list_phase_mfiles(stimdir, phaseName);
if isempty(mfiles)
    warning('[%s] No mfiles found in stim/', phaseName);
    return
end

t0_run = RUN.run.(phaseName).t0_ts;
t1_run = RUN.run.(phaseName).t1_ts;

trigRun = subset_trigOE_to_window(trigOE, t0_run, t1_run);
trigPhase = trigRun;

t = get_baphy_ttl_seconds(trigRun);
if isempty(t)
    warning('[%s] No trigRun.baphy.t_raw_s', phaseName);
    return
end

t = RA_prune_close_pulses(t, 0.5);

blk = RAExp_split_blocks(t, blockGap_s);
nBlk = size(blk,1);
nMf  = numel(mfiles);
nUse = min(nBlk, nMf);

if nUse < nMf || nUse < nBlk
    warning('[%s] blocks=%d, mfiles=%d. Using %d pairs in order.', phaseName, nBlk, nMf, nUse);
end

Bparts = cell(nUse,1);
Eparts = cell(nUse,1);

for k = 1:nUse
    mf = mfiles{k};
    infoMf = RA_get_mfile_info(mf);

    blockT = t(blk(k,1):blk(k,2));
    selT = RA_select_block_subset(blockT, infoMf.n_trials, infoMf.trial_span_s);

    t0_s = selT(1) - pad_s;
    t1_s = selT(end) + pad_s;

    trK = subset_trigOE_to_window(trigOE, round(t0_s*TsRate), round(t1_s*TsRate));
    trK = RA_overwrite_baphy_ttls(trK, selT);

    outTag = phaseName;
    if k > 1
        outTag = sprintf('%s_%d', phaseName, k-1);
    end

    Bk = RA_parse_baphy_active(datapath, trK, ...
        'MFile', mf, 'DLCSubdir', dlcSubdir, 'OutTag', outTag);

    Ek = RA_build_epochs_active(datapath, trK, Bk, dlcSubdir);

    parts(k).mfile = mf;
    parts(k).t0_s  = t0_s;
    parts(k).t1_s  = t1_s;
    parts(k).B     = Bk;
    parts(k).E     = Ek;
    parts(k).n_trials_expected = infoMf.n_trials;
    parts(k).n_ttl_used = numel(selT);

    Bparts{k} = Bk;
    Eparts{k} = Ek;
end

if isempty(Bparts)
    warning('[%s] No stitched parts were created.', phaseName);
    return
end

Bphase = RAExp_stitch_baphy(Bparts);
Ephase = RAExp_stitch_epochs(Eparts);

end


function mfiles = RA_list_training_mfiles(stimdir)
% Sort training mfiles as ..._1.m, ..._2.m ... ; fallback to name order.
d = dir(fullfile(stimdir, '*.m'));
if isempty(d), mfiles = {}; return; end

k = nan(numel(d),1);
for i = 1:numel(d)
    tok = regexp(d(i).name, '_(\d+)\.m$', 'tokens', 'once');
    if ~isempty(tok)
        k(i) = str2double(tok{1});
    end
end

if all(isnan(k))
    [~, ord] = sort(lower({d.name}));
else
    k(isnan(k)) = 0;
    [~, ord] = sort(k);
end
d = d(ord);

mfiles = cell(numel(d),1);
for i = 1:numel(d)
    mfiles{i} = fullfile(stimdir, d(i).name);
end
end

function t = get_baphy_ttl_seconds(trigS)
t = [];
if isfield(trigS,'baphy') && isstruct(trigS.baphy)
    if isfield(trigS.baphy,'t_raw_s') && ~isempty(trigS.baphy.t_raw_s)
        t = trigS.baphy.t_raw_s(:);
    elseif isfield(trigS.baphy,'trial_start_ts') && ~isempty(trigS.baphy.trial_start_ts)
        t = trigS.baphy.trial_start_ts(:) ./ 1e4;
    elseif isfield(trigS.baphy,'t_raw_ts') && ~isempty(trigS.baphy.t_raw_ts)
        t = trigS.baphy.t_raw_ts(:) ./ 1e4;
    end
end
end

function t2 = RA_prune_close_pulses(t, minISI_s)
t2 = t(:);
if numel(t2) < 2
    return
end
keep = [true; diff(t2) > minISI_s];
t2 = t2(keep);
end

function infoMf = RA_get_mfile_info(mfile)
infoMf = struct('n_trials', NaN, 'trial_span_s', NaN);
try
    clear globalparams exptparams exptevents
    run(mfile);
    if exist('exptevents','var')
        trAll = [exptevents.Trial];
        trAll = trAll(trAll > 0 & isfinite(trAll));
        if ~isempty(trAll)
            infoMf.n_trials = double(max(trAll));
        end

        if isfinite(infoMf.n_trials) && infoMf.n_trials > 0
            trId = [exptevents.Trial];
            stopMax = nan(infoMf.n_trials,1);
            for tt = 1:infoMf.n_trials
                kk = find(trId == tt);
                if isempty(kk), continue; end
                st = [exptevents(kk).StopTime];
                if ~isempty(st)
                    stopMax(tt) = max(st);
                end
            end
            if any(isfinite(stopMax))
                infoMf.trial_span_s = sum(stopMax(isfinite(stopMax)));
            end
        end
    end
catch
end
end

function selT = RA_select_block_subset(blockT, nExpected, trialSpan_s)
selT = blockT(:);
if isempty(selT)
    return
end
if ~isfinite(nExpected) || nExpected <= 0
    return
end
if numel(selT) <= nExpected
    selT = selT(1:min(numel(selT), nExpected));
    return
end

nWin = numel(selT) - nExpected + 1;
bestCost = inf;
bestIdx = 1;

for i = 1:nWin
    sub = selT(i:i+nExpected-1);
    cost = 0;

    dt = diff(sub);
    if ~isempty(dt)
        meddt = median(dt);
        if isfinite(meddt) && meddt > 0
            cost = cost + sum(dt < 0.25 * meddt);
        end
    end

    if isfinite(trialSpan_s) && trialSpan_s > 0
        spanErr = abs((sub(end) - sub(1)) - trialSpan_s) / trialSpan_s;
        cost = cost + spanErr;
    else
        cost = cost + i / max(1, nWin);
    end

    if cost < bestCost
        bestCost = cost;
        bestIdx = i;
    end
end

selT = selT(bestIdx:bestIdx+nExpected-1);
end

function trK = RA_overwrite_baphy_ttls(trK, selT_s)
selT_s = selT_s(:);
if ~isfield(trK,'baphy') || ~isstruct(trK.baphy)
    trK.baphy = struct();
end
trK.baphy.t_raw_s = selT_s;
trK.baphy.t_raw_ts = round(selT_s * 1e4);
trK.baphy.trial_start_ts = trK.baphy.t_raw_ts;
trK.baphy.n_trials = numel(selT_s);
end

function mfiles = RAExp_list_phase_mfiles(stimdir, phaseName)
% Sort: ..._<Phase>.m, ..._<Phase>_1.m, ..._<Phase>_2.m ...
d = dir(fullfile(stimdir, ['*_' phaseName '*.m']));
if isempty(d), mfiles = {}; return; end

k = zeros(numel(d),1);
for i = 1:numel(d)
    nm = d(i).name;
    tok = regexp(nm, ['_' phaseName '_(\d+)\.m$'], 'tokens','once');
    if isempty(tok)
        k(i) = 0;
    else
        k(i) = str2double(tok{1});
    end
end
[~,ord] = sort(k);
d = d(ord);

mfiles = cell(numel(d),1);
for i = 1:numel(d)
    mfiles{i} = fullfile(stimdir, d(i).name);
end
end

function blk = RAExp_split_blocks(t_s, gap_s)
% t_s are trial pulses (RA: 1 pulse per trial). Split by gaps > gap_s.
if isempty(t_s), blk = zeros(0,2); return; end
dt = diff(t_s(:));
cut = find(dt > gap_s);
starts = [1; cut+1];
stops  = [cut; numel(t_s)];
blk = [starts stops];
end

function Bout = RAExp_stitch_baphy(Bparts)
% Concatenate B.trial fields; adjust B.idx fields by offsets; keep per-part in Bout.part.
Bout = Bparts{1};
Bout.part = Bparts(:);

n = 0;
for i = 1:numel(Bparts)
    if isfield(Bparts{i},'n_trials'), n = n + Bparts{i}.n_trials; end
end
Bout.n_trials = n;

% stitch trial fields (numeric + logical + cell + string)
if isfield(Bout,'trial') && isstruct(Bout.trial)
    fn = fieldnames(Bout.trial);
    for f = 1:numel(fn)
        name = fn{f};

        vals = cell(numel(Bparts),1);
        okNumLog = true;
        okCell   = true;
        okStr    = true;

        for i = 1:numel(Bparts)
            if ~isfield(Bparts{i},'trial') || ~isfield(Bparts{i}.trial,name)
                okNumLog = false; okCell = false; okStr = false;
                break
            end

            v = Bparts{i}.trial.(name);

            okNumLog = okNumLog && (isnumeric(v) || islogical(v));
            okCell   = okCell   && iscell(v);
            okStr    = okStr    && isstring(v);

            % store as column for consistent vertcat
            if isnumeric(v) || islogical(v) || isstring(v)
                vals{i} = v(:);
            elseif iscell(v)
                vals{i} = v(:);
            else
                okNumLog = false; okCell = false; okStr = false;
            end
        end

        try
            if okNumLog
                Bout.trial.(name) = vertcat(vals{:});
            elseif okCell
                Bout.trial.(name) = vertcat(vals{:});
            elseif okStr
                Bout.trial.(name) = vertcat(vals{:});
            end
        catch
        end
    end
end

% ensure trial_id matches stitched length
Bout.trial_id = (1:Bout.n_trials)';

% ---------------- stitch Perf fields (per-trial vectors) ----------------
if isfield(Bout,'Perf') && isstruct(Bout.Perf)
    nParts = numel(Bparts);
    nPer = zeros(nParts,1);
    for i = 1:nParts
        nPer(i) = get_n_trials_safe(Bparts{i});
    end

    Pout = struct();
    Pout.n_trials = Bout.n_trials;

    pfn = fieldnames(Bparts{1}.Perf);
    for f = 1:numel(pfn)
        nm = pfn{f};

        vals = cell(nParts,1);
        okNumLog = true;
        okCell   = true;
        okStr    = true;
        okLen    = true;

        for i = 1:nParts
            if ~isfield(Bparts{i},'Perf') || ~isfield(Bparts{i}.Perf,nm)
                okNumLog = false; okCell = false; okStr = false; okLen = false;
                break
            end
            v = Bparts{i}.Perf.(nm);

            okNumLog = okNumLog && (isnumeric(v) || islogical(v));
            okCell   = okCell   && iscell(v);
            okStr    = okStr    && isstring(v);

            if (isnumeric(v) || islogical(v) || isstring(v) || iscell(v)) && isvector(v)
                okLen = okLen && (numel(v) == nPer(i));
                vals{i} = v(:);
            else
                okLen = false;
            end
        end

        % Only stitch fields that are per-trial vectors in *every* part.
        if okLen && (okNumLog || okCell || okStr)
            try
                Pout.(nm) = vertcat(vals{:});
            catch
            end
        end
    end

    % keep original per-part Perf for debugging (optional)
    Pout.part = cell(nParts,1);
    for i = 1:nParts
        if isfield(Bparts{i},'Perf'), Pout.part{i} = Bparts{i}.Perf; else, Pout.part{i} = []; end
    end

    Bout.Perf = Pout;
end

% stitch idx (shift ONLY numeric trial-index vectors; concatenate others)
if isfield(Bout,'idx') && isstruct(Bout.idx)
    idxNames = fieldnames(Bout.idx);
    for j = 1:numel(idxNames)
        nm = idxNames{j};

        accNum  = [];
        accLog  = [];
        accCell = {};
        accStr  = [];

        off = 0;
        for i = 1:numel(Bparts)
            if ~isfield(Bparts{i},'idx') || ~isfield(Bparts{i}.idx,nm) || isempty(Bparts{i}.idx.(nm))
                off = off + get_n_trials_safe(Bparts{i});
                continue
            end

            v = Bparts{i}.idx.(nm);

            if isnumeric(v)
                tmp = v(:) + off;                 % shift numeric indices
                accNum = [accNum; tmp]; %#ok<AGROW>

            elseif islogical(v)
                accLog = [accLog; v(:)]; %#ok<AGROW> % concatenate logical mask

            elseif iscell(v)
                accCell = [accCell; v(:)]; %#ok<AGROW> % concatenate labels/cells (NO shift)

            elseif isstring(v)
                accStr = [accStr; v(:)]; %#ok<AGROW>   % concatenate strings (NO shift)

            else
                % fallback: try treat as cell
                try
                    vv = cellstr(v);
                    accCell = [accCell; vv(:)]; %#ok<AGROW>
                catch
                end
            end

            off = off + get_n_trials_safe(Bparts{i});
        end

        if ~isempty(accNum)
            Bout.idx.(nm) = accNum;
        elseif ~isempty(accLog)
            Bout.idx.(nm) = accLog;
        elseif ~isempty(accStr)
            Bout.idx.(nm) = accStr;
        elseif ~isempty(accCell)
            Bout.idx.(nm) = accCell;
        else
            Bout.idx.(nm) = [];
        end
    end
end

end 

function n = get_n_trials_safe(B)
n = 0;
if isfield(B,'n_trials') && ~isempty(B.n_trials) && isnumeric(B.n_trials)
    n = B.n_trials;
elseif isfield(B,'trial') && isstruct(B.trial)
    fn = fieldnames(B.trial);
    for k = 1:numel(fn)
        v = B.trial.(fn{k});
        if isnumeric(v) || iscell(v) || isstring(v)
            n = numel(v);
            break
        end
    end
end
end
function Eout = RAExp_stitch_epochs(Eparts)
% Union intervalSets and concatenate vector fields.
Eout = Eparts{1};
Eout.part = Eparts(:);

fn = fieldnames(Eout);
for f = 1:numel(fn)
    name = fn{f};
    if strcmp(name,'part'), continue; end
    
    v0 = Eout.(name);
    
    if isa(v0,'intervalSet')
        u = v0;
        for i = 2:numel(Eparts)
            if isfield(Eparts{i},name) && isa(Eparts{i}.(name),'intervalSet')
                u = union(u, Eparts{i}.(name));
            end
        end
        Eout.(name) = u;
        
    elseif isnumeric(v0) && isvector(v0)
        vals = cell(numel(Eparts),1);
        ok = true;
        for i = 1:numel(Eparts)
            if ~isfield(Eparts{i},name) || ~isnumeric(Eparts{i}.(name))
                ok = false; break
            end
            vals{i} = Eparts{i}.(name)(:);
        end
        if ok
            Eout.(name) = vertcat(vals{:});
        end
    end
end
end

function n = get_baphy_n_trials_any(B)
n = [];
if isstruct(B) && isfield(B,'n_trials') && isnumeric(B.n_trials) && ~isempty(B.n_trials)
    n = B.n_trials;
    return
end
if isstruct(B) && isfield(B,'Conditioning') && isfield(B.Conditioning,'n_trials')
    n = 0;
    n = n + double(B.Conditioning.n_trials);
    if isfield(B,'PostTest') && isfield(B.PostTest,'n_trials')
        n = n + double(B.PostTest.n_trials);
    end
    return
end
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