function React_Passive_Behaviour_preprocessing()

% This is a master script to process the DLC behavioural ferret data for
% ReactActive project
% Steps:
%   - Synchronization with LFP signal (sync_behaviour_ephys). It creates a correct timeline taking into account the delay between the video and ephys
%   - Generation of the basic figures (behaviour_analysis)
% Under construction  - Study the correlation between OB/Cortical/Hippocampal gamma and pupil area (gamma_pupil_corr)
% Under construction  - Producing the composition video with all variables synced (composition_video_OB_DLC_ferret)

k = 1;
session_dlc = {};

for c = 1:length(sessions)
    dlc_path = fullfile(sessions{c}, 'video'); 
    files = dir(fullfile(dlc_path, '*_filtered.csv')); % Search for files ending with "_filtered.csv"
    
    if ~isempty(files) % Check if there are any matching files
        session_dlc{k} = sessions{c}; % Store the session path
        k = k + 1;
    else
        disp([Dir{selection}.path{c} ' - No DLC found']);
    end
end

session_dlc = session_dlc';

%% -------------------------------------- PREPROCESSING -------------------------------------- 
% SESS: Do the basic DLC pre-processing
for sess = 1:length(session_dlc)
    disp('Analysing DLC data...')
    disp([num2str(length(session_dlc)-sess + 1) '/' num2str(length(session_dlc)) ' left'])
    disp(['Running session: ' session_dlc{sess}])
    behaviour_preprocessing(session_dlc{sess})
end

% SESS: Check the quality of tracking on a short episode
% range = [8 50];
% marker = {'pupil_area_007', 'cheek_center_mvt', 'nostril_center_mvt', 'jaw_center_mvt', 'tongue_center_mvt', 'spout_likelihood', 'jaw_center'};
% for sess = 1:length(session_dlc)
%     disp(['Running session: ' session_dlc{sess}])
%     disp([num2str(length(session_dlc)-sess + 1) '/' num2str(length(session_dlc)) ' left'])
%     for m = 2:numel(marker)
%         disp(['Running marker: ' marker{m}])
%         disp([num2str(numel(marker)-m + 1) '/' num2str(numel(marker)) ' left'])
%         tracking_check(session_dlc{sess}, marker{m}, range) %generates movie to verify the sync and accuracy of tracked markers and video
%     end
% end

% SESS: Extract the behavioural data & Create trial epochs
for sess = 1:length(session_dlc)
    disp('Analysing Baphy data...')
    disp([num2str(length(session_dlc)-sess + 1) '/' num2str(length(session_dlc)) ' left'])    
    disp(['Running session: ' session_dlc{sess}])
    baphy_preprocessing(session_dlc{sess})
end
end