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

    % Synchronize LFP and DLC ; Produces synced timeline in DLC_data.mat
    if opts.do_sync
        disp('Syncing DLC and Ephys...')
        sync_behaviour_ephys(datapath);
    end

    % Do the basic DLC pre-processing
    if opts.do_basic_preproc
        disp('Basic DLC preprocessing...')
        behaviour_preprocessing(datapath);
    end
    
    % Do motion energy analysis
    if opts.do_motion_energy
        motion_energy(session_dlc{sess})
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