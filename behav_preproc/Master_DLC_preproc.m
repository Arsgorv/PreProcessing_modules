function Master_DLC_preproc(sessions, opts)

% This is a master script to preprocess the DLC behavioural ferret data
% Steps:
%   - Synchronization with LFP signal (sync_behaviour_ephys). It creates a correct timeline taking into account the delay between the video and ephys
%   - Do the basic DLC pre-processing
%   - Calculate motion energy
%   - Check the quality of DLC tracking

%% session-based DLC preprocessing
if nargin < 2, opts = struct(); end
if ~isfield(opts,'do_sync'), opts.do_sync = true; end
if ~isfield(opts,'do_basic_preproc'), opts.do_basic_preproc = true; end
if ~isfield(opts,'do_motion_energy'), opts.do_motion_energy = false; end

session_dlc = filter_sessions_with_dlc(sessions);

for sess = 1:length(session_dlc)
    datapath = session_dlc{sess};
    disp(['[behav_preproc] ' datapath])
    disp([num2str(numel(session_dlc)-sess+1) '/' num2str(numel(session_dlc)) ' left'])
    
    isRAexp = ~isempty(dir(fullfile(datapath,'ephys','*_RA_PreSleep*'))) || ...
        ~isempty(dir(fullfile(datapath,'ephys','*_RA_PostTest*')));
    
    if ~isRAexp
        if opts.do_sync
            disp('Syncing DLC and Ephys...')
            sync_behaviour_ephys(datapath);
        end
        if opts.do_basic_preproc
            disp('Basic DLC preprocessing...')
            behaviour_preprocessing(datapath); 
        end
        if opts.do_motion_energy
            disp('Do motion energy preprocessing...')
            motion_energy(datapath)
        end
    else
        if opts.do_sync
            % Synchronize LFP and DLC ; Produces synced timeline in DLC_data.mat
            disp('Syncing DLC and Ephys...')
            sync_behaviour_ephys(datapath,'Conditioning');
            sync_behaviour_ephys(datapath,'PostTest');
        end
        if opts.do_basic_preproc
            % Do the basic DLC pre-processing
            disp('Basic DLC preprocessing...')
            behaviour_preprocessing(datapath,'Conditioning');
            behaviour_preprocessing(datapath,'PostTest');
        end
        if opts.do_motion_energy
            % Do motion energy analysis            
            disp('Do motion energy preprocessing...')
            motion_energy(datapath,'Conditioning')
            motion_energy(datapath,'PostTest')
        end
    end

    
    % Check the quality of tracking on a short episode
    % range = [8 50];
    % marker = {'pupil_area_007', 'cheek_center_mvt', 'nostril_center_mvt', 'jaw_center_mvt', 'tongue_center_mvt', 'spout_likelihood', 'jaw_center'};
    %     for m = 2:numel(marker)
    %         disp(['Running marker: ' marker{m}])
    %         disp([num2str(numel(marker)-m + 1) '/' num2str(numel(marker)) ' left'])
    %         tracking_check(session_dlc{sess}, marker{m}, range) %generates movie to verify the sync and accuracy of tracked markers and video
    %     end
    
    close all
end

end