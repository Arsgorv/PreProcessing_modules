function run_PSTH_analysis_all_sess
sessions = {...
    '20241205_TORCs', ...
    '20241206_TORCs', ...
    '20241209_TORCs', ...
    '20241210_TORCs', ...
    '20241211_TORCs', ...
    '20241212_TORCs', ...
    '20241213_TORCs', ...
    '20241214_TORCs', ...
    '20241221_TORCs', ...
    '20241223_TORCs'};


exp_path  = '/media/nas8/OB_ferret_AG_BM/Shropshire/freely-moving/';
% session = sessions{5}; %5

for sess = 1:length(sessions)
    tic
    session = sessions{sess};
    directory = [exp_path session];
    cd(directory)
    disp(session)
    rasterDir = fullfile([directory '/wave_clus/'], 'raster_figures');
    if ~exist(rasterDir, 'dir'), mkdir(rasterDir); end
    
    
    %% SECTION 1.2: Load wave_clus Spikes (or KlustaKwik if you prefer)
    
    [spikes, metadata] = load_wave_clus(directory);
    
    % Optional: Remove bad data for a specific session
    % if strcmp(session, '20241205_TORCs')
    %     metadata([64:75], :) = [];
    %     spikes(:, [64:75])   = [];
    % end
    
    %% SECTION 1.3: Load LFP Data & Extract Stimulus Onsets
    lfpFile = fullfile(directory, 'LFPData', 'LFP111.mat');
    if exist(lfpFile,'file')
        load(lfpFile,'LFP');
    else
        error('LFP file not found: %s', lfpFile);
    end
    
    % Extract stimulus onsets:
    StartSound = thresholdIntervals(LFP, 4000, 'Direction','Above');
    Starttime  = Start(StartSound);
    
    % set(0, 'DefaultFigureWindowStyle','docked');  % or 'normal'
    
    % Convert to convinient format for Yves (STRF)
    % stim_start = Starttime/1e4;
    % save([exp_path 'stim_starts/' 'stim_start_' session], 'stim_start')
    
    %% WILL BE INTEGRATED WITH SECTION 1.4 SOON...... Selecting “Good” Clusters (colIndex) Based on Preselection
    clear chcl
    chcl = get_chcl(session);  % function that returns list of [channel cluster]
    colIndex = [];
    for s = 1:length(chcl)
        % Find row in metadata matching [channel cluster]
        row = find(metadata(:,1)==chcl{s}(1) & metadata(:,2)==chcl{s}(2));
        colIndex(s) = row;
    end
    % spikeTimes = spikes(:, colIndex);
    % spikeTimes = spikeTimes(~isnan(spikeTimes)); % Remove NaN values
    
    %% SECTION 1.5: Create a Cell Array A{n} for Neurons
    % A{n} = ts(...) for each selected cluster
    A = {};
    n=1;
    for i = colIndex
        A{n} = ts(spikes(:,i)*1e1);
        n = n+1;
    end
    
    %% SECTION 1.6: Load Sleep & Gamma Data
    % Load SleepScoring info, e.g. Epoch, Wake, Sleep, SWSEpoch, REMEpoch, etc.
    sleepFile = fullfile(directory,'SleepScoring_OBGamma.mat');
    if exist(sleepFile,'file')
        load(sleepFile,'Wake','Sleep','SWSEpoch','REMEpoch','Epoch_S1','Epoch_S2','SmoothGamma','Epoch');
    else
        warning('No SleepScoring_OBGamma.mat found in %s', directory);
    end
    
    %% SECTION 4.1: Build Stimulus Onset States for PSTH
    %     stim_onset = ts(Starttime);
    %     stim_onset_Wake = Restrict(stim_onset, Wake);
    %     stim_onset_Sleep = Restrict(stim_onset, Sleep);
    %     stim_onset_NREM = Restrict(stim_onset, SWSEpoch);
    %     stim_onset_REM = Restrict(stim_onset, REMEpoch);
    %
    %     stim_onset_NREM1 = Restrict(stim_onset, (Sleep & Epoch_S2) - REMEpoch);
    %     stim_onset_NREM2 = Restrict(stim_onset, (Sleep & Epoch_S1) - REMEpoch);
    %
    %     stim_states = {stim_onset_Wake; stim_onset_Sleep; stim_onset_NREM; ...
    %         stim_onset_REM; stim_onset_NREM1; stim_onset_NREM2};
    %     conditionNames = {'Wake','Sleep','NREM','REM','NREM1','NREM2'};
    
    %% SECTION 4.2: PSTH ANALYSIS (Raw & Z-score) + Population Figures + Controls
%         fig_vis = 'off';
    %     [rawRates, popRawPSTH, zRates, popZscoredPSTH, zBaselineCorrectedRates, popZscoredPSTH_BC, timeCenters] = PSTH_Analysis_AllNeurons(A, chcl, stim_states, conditionNames, directory, fig_vis);
    %
%         outFile = fullfile([directory '/wave_clus/'], sprintf('spike_analysis_%s.mat', session));
    %     save(outFile, ...
    %         'rawRates', 'zRates', 'zBaselineCorrectedRates',...
    %         'popRawPSTH', 'popZscoredPSTH', 'popZscoredPSTH_BC', ...
    %         'timeCenters', 'conditionNames', ...
    %         'metadata', 'colIndex', 'chcl', ...
    %         '-v7.3');
    %     fprintf('Saved PSTH results to %s\n', outFile);
    %
    %% FIGURE: WakeGamma-dependent PSTH
    fig_vis = 'off';
    
    [rawRatesGamma, popRawPSTHGamma, zRatesGamma, popZscoredPSTHGamma, zBaselineCorrectedRatesGamma, popZscoredPSTH_BCGamma, timeCentersGamma, GammaEpochs, GammaConditionNames] = PSTH_GammaBins_SingleNeuron(A, chcl, SmoothGamma, Wake, Starttime, directory, fig_vis);
    outFile = fullfile([directory '/wave_clus/'], sprintf('spike_analysis_%s.mat', session));
    
    save(outFile, ...
        'rawRatesGamma', 'zRatesGamma', 'zBaselineCorrectedRatesGamma',...
        'popRawPSTHGamma', 'popZscoredPSTHGamma', 'popZscoredPSTH_BCGamma', ...
        'timeCentersGamma', 'GammaEpochs', 'GammaConditionNames', ...
        'metadata', 'colIndex', 'chcl', ...
        '-v7.3', '-append');
    fprintf('Saved PSTH results to %s\n', outFile);
    
    toc/60
end
end