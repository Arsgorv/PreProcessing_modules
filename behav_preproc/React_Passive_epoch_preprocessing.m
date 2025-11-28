function React_Passive_epoch_preprocessing
% Epochs:
%{
1. Define common timeline (extract triggers ; extract stimulus identity)

2. Create Epochs
- fUS epoch
- Trial epoch
- Stimulus_1 epoch
- Stimulus_2 epoch

3. apply fUS denoising pipeline and define fUS noise epoch
%}



%% Extract triggers to sync fUS, baphy and OE
if contains(datapath, 'Edel') || contains(datapath, 'Chabichou')
    
    % Load trigger information
    
    list_m = dir([datapath filesep 'baphy' filesep '*.m']);
    m_file_name = list_m(1).name;
    run([datapath filesep 'baphy' filesep m_file_name]);
    
    list_csv = dir([datapath filesep 'baphy' filesep '*.csv']);
    trig_file_name = list_csv(1).name;
    trig_file = readtable([datapath filesep 'baphy' filesep trig_file_name]);
    
    filenames = [datapath filesep 'TrigFiles' filesep '*REA_' params.slice{sess} '_' params.pair{sess} '_' params.session{sess} '*.csv'];
    listings = dir(filenames);
    
    if isempty(listings)
        filenames = [datapath '/TrigFiles/*REA_' params.slice{sess} '_p' params.pair{sess} '_' params.session{sess} '*.csv'];
        listings = dir(filenames);
    end
    
    %         trig_file = listings(1).name;
    %         file = readtable([datapath filesep 'TrigFiles' filesep trig_file]);
    
    % Extract trigs and cut to trials
    % fUS trigs
    trig_info = table2cell(trig_file(:,2));
    f_trigs = trig_info(cellfun(@(x) isequal(x,1),strfind(trig_info,'F-')));
    f_trigs = cell2mat(cellfun(@(x) str2double(x(3:end)),f_trigs,'UniformOutput',false));
    
    % baphy trigs
    b_trigs = trig_info(cellfun(@(x) isequal(x,1),strfind(trig_info,'B-')));
    b_trigs = cell2mat(cellfun(@(x) str2double(x(3:end)),b_trigs,'UniformOutput',false));
    
    % fix issue when first fUS frame had a non-0 timing (not supposed to happen but is there sometimes)
    b_trigs = b_trigs - f_trigs(1);
    f_trigs = f_trigs - f_trigs(1);
    
elseif contains(datapath, 'Kosichka')
    
    % SESS: Synchronize LFP and DLC ; Produces synced timeline in DLC_data.mat
    for sess = 8:length(session_dlc)
        disp(['Running session: ' session_dlc{sess}])
        disp([num2str(length(session_dlc)-sess + 1) '/' num2str(length(session_dlc)) ' left'])
        disp('Syncing DLC and Ephys...')
        sync_behaviour_ephys(session_dlc{sess})
    end
    
end

%% Cut data_cat into trials
% Cut into trials
cut_in_trials_all_sess

% easy way to check that I have irregular fUS sampling
% figure;plot(diff(f_trigs))

%         cut_into_trials_FL
%         cut_into_trials_AG


% cut into trials
[rawdatacut, trial_timings] = cut_into_trials_AB(data_cat, sess, f_trigs, b_trigs, plt, exp_info);
%         [rawdatacut, trial_timings] = cut_into_trials_AB(raw_data, sess, f_trigs, b_trigs, plt, exp_info);

data_cut_in_trials(:, :, :, :, sess) = rawdatacut;
n_trials = size(rawdatacut, 3);
save('data_cut_in_trials_PreExp', 'data_cut_in_trials', 'trial_timings');








end