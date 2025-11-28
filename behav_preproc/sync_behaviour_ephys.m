function sync_behaviour_ephys(datapath)
%%%%%% Synchronization of video and OB (LFP) signal %%%%%%%%%
% This function:
%   1. Loads LFP data and extracts trigger timestamps.
%   2. Determines the time of the first trigger (time_1st_trig).
%   3. Loads video timestamps (data_csv) and DLC tracking data.
%   4. Checks if DLC frame count equals video timestamps; if not, interpolates DLC.
%   5. Computes delay between openephys and video start and aligns DLC data.

%% ----------- Part 1: Load and Process LFP Data (Openephys) -----------
% Load LFP data (here using different paths based on your datapath)
Session_params.plt = 0;
if contains(datapath, 'Shropshire')
    l = 112;
    load(fullfile(datapath, 'LFPData', ['LFP' num2str(l) '.mat']));
    LFP_ferret{l} = LFP;
elseif contains(datapath, 'Brynza')
    l = 42;
    load(fullfile(datapath, 'LFPData', ['LFP' num2str(l) '.mat']));
    LFP_ferret{l} = LFP;
elseif contains(datapath, 'Labneh')
    l = 42;
    load(fullfile(datapath, 'LFPData', ['LFP' num2str(l) '.mat']));
    LFP_ferret{l} = LFP;
elseif contains(datapath, 'Tvorozhok')
    face_04_lfp = 23; 
    eye_07_lfp = 24; 
    load(fullfile(datapath, 'LFPData', ['LFP' num2str(face_04_lfp) '.mat']));
    LFP_ferret{face_04_lfp} = LFP;
    load(fullfile(datapath, 'LFPData', ['LFP' num2str(eye_07_lfp) '.mat']));
    LFP_ferret{eye_07_lfp} = LFP;
elseif contains(datapath, 'Kosichka')
    face_04_lfp = 23;
    eye_07_lfp = 24; 
    load(fullfile(datapath, 'LFPData', ['LFP' num2str(face_04_lfp) '.mat']));
    LFP_ferret{face_04_lfp} = LFP;
    load(fullfile(datapath, 'LFPData', ['LFP' num2str(eye_07_lfp) '.mat']));
    LFP_ferret{eye_07_lfp} = LFP;
end

% Define the time scales
time_lfp = Range(LFP_ferret{face_04_lfp}); % in ts
data_lfp = {Data(LFP_ferret{face_04_lfp}) ; Data(LFP_ferret{eye_07_lfp})};

for idx = 1:2
    % Process LFP trigger channel to extract trigger timestamps.
    LFP_trig_data = data_lfp{idx};
    tmp = diff(LFP_trig_data);
    threshold_diff = max(tmp)/3;  % you might adjust this threshold as needed
    rise_fall_idx = find(tmp >= threshold_diff);
    rise_fall = tmp(rise_fall_idx);
    
    % Set a minimum gap (in indices) to separate events:
    threshold_event = 4;
    num_points = numel(rise_fall_idx);
    
    % Determine event indices (start of a new event)
    event_indices = 1;
    for i = 2:num_points
        if (rise_fall_idx(i) - rise_fall_idx(i-1)) > threshold_event
            event_indices = [event_indices, i];
        end
    end
    event_indices = [event_indices, num_points];
    
    % For each event, pick the peak value and its index.
    num_events = numel(event_indices) - 1;
    peak_values = nan(num_events, 1);
    peak_indices = nan(num_events, 1);
    for i = 1:num_events
        event_data = LFP_trig_data(rise_fall_idx(event_indices(i):event_indices(i+1)-1));
        if ~isempty(event_data)
            [~, max_idx] = max(event_data);
            current_peak_index = rise_fall_idx(event_indices(i) - 1 + max_idx);
            peak_values(i) = event_data(max_idx);
            peak_indices(i) = current_peak_index;
        end
    end
    
    % for the session where I had an outlier as the first index
    % peak_values = nan(num_events-1, 1);
    % peak_indices = nan(num_events-1, 1);
    
    % (Handle possible NaN in last event)
    if any(isnan(peak_indices))
        idxNaN = find(isnan(peak_indices));
        peak_indices(idxNaN) = peak_indices(end-1) + round(nanmean(diff(peak_indices)));
        peak_values(idxNaN) = peak_values(end-1) + round(nanmean(diff(peak_values)));
        disp('Last value(s) of peak_indices/peak_values were NaN. Replaced with approximated values.');
    end
    
    % Remove outlier trigger(s) if needed.
    % Inter-peak interval and outliers (modify the threshold)
    % Calculate interpeak intervals
    interpeak_intervals = diff(peak_indices);
    % Calculate z-scores for the interpeak intervals
    z_scores = zscore(interpeak_intervals);
    % Set the threshold for outlier detection (e.g., z-score > 3 or < -3)
    threshold = 3.5;
    % Find the outlier intervals
    outlier_indices = find(z_scores > threshold | z_scores < -threshold);
    outlier_intervals = interpeak_intervals(outlier_indices);
    if ismember(1, outlier_indices)
        % Remove outliers from peak_indices and peak_values
        peak_indices(outlier_indices == 1) = [];
        peak_values(outlier_indices == 1) = [];
        %
        %     peak_indices(outlier_indices+1) = [];
        %     peak_values(outlier_indices) = [];
        outlier_indices(1) = [];
    end
    
    if outlier_indices == size(peak_indices, 1)
        peak_indices(end) = [];
        peak_values(end) = [];
    end
    
    % Define the first trigger time from LFP as the baseline.
    time_1st_trig(idx) = time_lfp(peak_indices(1));
    if idx == 1
        fprintf('FACE cam: first openephys trigger time: %fs (in LFP time scale)\n', time_1st_trig(idx)/1e4);
    else
        fprintf('EYE cam: first openephys trigger time: %fs (in LFP time scale)\n', time_1st_trig(idx)/1e4);
    end
    % For synchronization, use the trigger timestamps.
    time_trig{idx} = time_lfp(peak_indices);  % these are your openephys trigger times in ts
 
    
    %% ----------- Part 2: Load DLC and Video Timestamp Data -----------
    % Set the DLC folder path
    dlc_path = fullfile(datapath, 'video');
    
    % Load video timestamps from the frames CSV (ground truth)
%     video_file = dir(fullfile(dlc_path, '*frames.csv'));
%     if isempty(video_file)
%         error('No video frames CSV file found in %s', dlc_path);
%     end
%     video_filename = video_file(1).name;
%     disp(['Video frames file: ' video_filename]);
%     data_csv = csvread(fullfile(dlc_path, video_filename));
%     
%     nCSV = length(data_csv);
%     fprintf('Number of video timestamps: %d\n', nCSV);
    
    % Load DLC tracking data (CSV file produced by DLC)
    dlc_file = dir(fullfile(dlc_path, '*_filtered.csv')); % Using the filtered (smoothed) file
    if isempty(dlc_file)
        error('No DLC filtered CSV file found in %s', dlc_path);
    end
    if contains(dlc_file(idx).name, '24934004')
        dlc_filename = dlc_file(1).name;
    elseif contains(dlc_file(idx).name, '24934007')
        dlc_filename = dlc_file(idx).name;
    end
    
    disp(['DLC data file: ' dlc_filename]);
    data = csvread(fullfile(dlc_path, dlc_filename), 3);
    [nFrames, nDims] = size(data);
    fprintf('DLC data has %d frames and %d dimensions.\n', nFrames, nDims);
    disp(['# OpenEphys triggers - # DLC frames = ' num2str(length(time_trig{idx})-nFrames) '. If it is a positive number, we will interpolate DLC data in the next step'])
    
%     if nCSV < nFrames
%         fprintf('DLC data is %d frame longer than the .csv. Trying to fix it...\n', nFrames-nCSV);
%         if data(end, 2) == 0
%             disp('The last DLC coordinate is 0. Removing it...')
%             data(end, :) = [];
%             [nFrames, nDims] = size(data);
%         end
%     end

    %% ----------- Part 3: Interpolate DLC Data (if needed) -----------
    % Our goal: ensure that DLC tracking data has one row per video timestamp.
    % There is no way to know when exactly camera droped the frames, so we make
    % an important assumption that they are dropped uniformly along the whole
    % recording. So, we evenly distribute the interpolation, which may be
    % completely wrong.
    nCSV = length(time_trig{idx});
    if nCSV ~= nFrames
        fprintf('Mismatch detected: OpenEphys triggers = %d, DLC frames = %d.\n', nCSV, nFrames);
        % Create a mapping from DLC frames to video timestamps.
        savedIndices = round(linspace(1, nCSV, nFrames));
        
        % (Optional) Check for duplicate indices
        if numel(unique(savedIndices)) < numel(savedIndices)
            warning('Duplicate indices detected in mapping. Check the linspace mapping.');
        end
        
        % Preallocate an array for the interpolated DLC data.
        dlcDataInterp = nan(nCSV, nDims);
        for col = 1:nDims
            % Interpolate each column using the openephys triggers timestamps corresponding to saved DLC frames.
            dlcDataInterp(:, col) = interp1(time_trig{idx}(savedIndices), data(:, col), time_trig{idx}, 'linear', 'extrap');
        end
        data = dlcDataInterp;
        fprintf('DLC data interpolated to %d entries (matching video timestamps).\n', size(data,1));
    else
        fprintf('Frame counts match; no interpolation needed.\n');
        savedIndices = (1:nCSV)';
    end
    
    %% ----------- Part 4: Synchronize DLC Data with OpenEphys Time -----------
    % Since openEphys (LFP) recording starts before video, we can compute a delay.
    % Use the first video timestamp (data_csv(1)) and the first LFP trigger (time_1st_trig).
%     delay = data_csv(1) - time_1st_trig(idx)/1e4;
%     fprintf('Computed delay between video and openephys: %f s\n', delay);
    
    % Align the video timestamps to openephys time:
%     videoTimestampsAligned = data_csv - delay;
%     time = ts(sort(videoTimestampsAligned)*1e4);
    
    % Combine the synchronized time with DLC tracking data.
    synchronizedData = [time_trig{idx}, data];
    if contains(dlc_file(idx).name, '24934004')
        csvwrite(fullfile(dlc_path, 'synchronized_DLC_data_face.csv'), synchronizedData);
    elseif contains(dlc_file(idx).name, '24934007')
        csvwrite(fullfile(dlc_path, 'synchronized_DLC_data_eye.csv'), synchronizedData);
    end
    % Save the synchronized DLC data for further analysis.
    
    fprintf('Synchronized DLC data saved to %s\n', fullfile(dlc_path, 'synchronized_DLC_data.csv'));
    
    %% ----------- Part 5: Additional Checks on Trigger vs. DLC Length -----------
    % Here we compare the number of openephys triggers (time_trig) with the video/DLC data length.
    if length(time_trig{idx}) > size(data,1)
        time_trig{idx} = time_trig{idx}(1:size(data,1));
        disp(['OpenEphys trigger signal is still longer than DLC data by ' num2str(length(time_trig{idx}) - size(data,1)) ' frames.']);
    elseif length(time_trig{idx}) < size(data,1)
        disp(['OpenEphys trigger signal is shorter than DLC data by ' num2str(size(data,1) - length(time_trig{idx})) ' frames. REVIEW THIS SESSION']);
    end
    fprintf('\n \n \n')        
    
    %% OLD: Threshold intervals version (doesn't work well) in ts
    % time_lfp = Range(LFP_ferret{l});
    %
    % % trig_onset = thresholdIntervals(LFP_ferret{l}, 1e4, 'Direction','Above');
    %
    % trig_onset = Start(thresholdIntervals(LFP_ferret{l}, 1e4, 'Direction','Above'));
    % if trig_onset(1) ~= 0
    %     time_1st_trig = trig_onset(1);
    % elseif trig_onset(1) == 0
    %     trig_onset(1) = [];
    %     time_1st_trig = trig_onset(1);
    % end
    %
    % disp([num2str(time_1st_trig/1e4) ' s'])
    %
    % [found, indices] = ismember(trig_onset, time_lfp);
    % peak_indices_new = indices(found);
    %
    % time_trig = trig_onset;
    %
    % % time_1st_trig_new = time_lfp(peak_indices_new(1));
    % % save([datapath '/DLC/DLC_data.mat'], 'time_1st_trig', 'time_trig', '-append');
    %
    % %YB
    % % time_trig = time_lfp(peak_indices)-time_1st_trig;
    % % time_trig = time_lfp(peak_indices);
    %
    % %
    % % % Save data
    % if exist([datapath '/DLC/DLC_data.mat'])
    %     save([datapath '/DLC/DLC_data.mat'], 'time_1st_trig', 'time_trig', '-append');
    % else
    %     save([datapath '/DLC/DLC_data.mat'], 'time_1st_trig', 'time_trig');
    % end
    %
    % missing_values  = [];
    % for i = 1:length(B)
    %     if all(abs(A-B(i))> tolerance)
    %         missing_values = [missing_values; i];
    %     end
    % end
    
    %% OLD: Control plots:
    
    %% ----------- (Optional) Control Plots -----------
    
    % Control plots
    if Session_params.plt(1) == 1
        % Plot the first tracking coordinate to check interpolation quality.
        figure;
        videoTimestampsAligned = time_trig{idx};
        plot(videoTimestampsAligned, data(:,2), 'b-', 'LineWidth', 1.5); hold on;
        % Overplot the original saved DLC data points (mapped via savedIndices)
        plot(videoTimestampsAligned(savedIndices), data(savedIndices,2), 'ro', 'MarkerSize', 6);
        xlabel('Time (s) [Aligned to openephys]');
        ylabel('Tracking Coordinate (Dimension 1)');
        legend('Interpolated DLC Data', 'Original Saved DLC Data');
        title('DLC Tracking Data Interpolation and Synchronization');
        grid on;
        
        %     % Plot the triggers with peak indices
        %     figure;
        %     subplot(2, 1, 1);
        %     plot(time_lfp, LFP_trig_data);
        %     hold on;
        %     scatter(time_lfp(rise_fall_idx), LFP_trig_data(rise_fall_idx), 'r', 'filled');
        %     xlabel('Index');
        %     ylabel('Value');
        %     title('Triggers with Rise-Fall Indexes');
        %
        %     subplot(2, 1, 2);
        %     plot(time_lfp, LFP_trig_data);
        %     hold on;
        %     scatter(time_lfp(peak_indices), peak_values, 'r', 'filled');
        %     xlabel('Index');
        %     ylabel('Value');
        %     title('Triggers with Peak Values');
        %
        %     % Display the plot
        %     grid on;
        
        % Plot the histogram of interpeak intervals with outliers
        figure;
        histogram(interpeak_intervals, 10000);
        hold on;
        histogram(outlier_intervals);
        legend('Interpeak Intervals', 'Outlier Intervals');
        xlabel('Interpeak Intervals');
        ylabel('Frequency');
        title('Histogram of Interpeak Intervals with Outliers');
        grid on;
        
        % Check the outliers and remove them (needs to be modified)
        %     % Plot the triggers with peak indices and outliers
        %     figure;
        %     plot(time_lfp, LFP_trig_data);
        %     hold on;
        %     scatter(rise_fall_idx, LFP_trig_data(rise_fall_idx), 'r', 'filled');
        %     scatter(time_lfp(peak_indices), peak_values, 'r', 'filled');
        %     scatter(time_lfp(peak_indices(outlier_indices)), peak_values(outlier_indices), 'g', 'filled');
        %     xlabel('Index');
        %     ylabel('Value');
        %     title('Triggers with Peak Indices (Outliers in Green)');
        %     grid on;
        
        figure; plot(LFP_trig_data); hold on; scatter(peak_indices, peak_values, 'r', 'filled')
        
        figure;
        plot(time_lfp, LFP_trig_data);
        grid on;
        hold on;
        % scatter(rise_fall_idx, LFP_trig_data(rise_fall_idx), 'r', 'filled');
        scatter(time_lfp(peak_indices), peak_values, 'r', 'filled');
        scatter(time_lfp(peak_indices(outlier_indices)), peak_values(outlier_indices), 'g', 'filled');
        scatter(time_trig{idx}, peak_values, 'b', 'filled');
        
        xlabel('Time');
        ylabel('Value');
        title('Mismatch between the csv timestamps and the timestamps, extracted from the trigger''s LFP idx');
        % Display the plot
        grid on;
        
        % sync DLC and triggers
        %     figure;
        %     grid on;
        %     subplot(4, 1, 1)
        %     plot(time_lfp, Data(LFP_OB))
        %     title('OB 24th channel')
        %
        %     subplot(4, 1, 2)
        %     plot(time_lfp, LFP_trig_data)
        %     hold on
        %     scatter(time_trig + time_1st_trig, peak_values, 'r', 'filled')
        %     title('LFP video triggers')
        %
        %     subplot(4, 1, 3)
        %     plot(time_trig + time_1st_trig, pupil_center(:, 1));
        %     hold on
        %     plot(time_trig + time_1st_trig, pupil_center(:, 1), '.', 'MarkerSize',15);
        %     title('DLC pupil center_x')
        %
        %     subplot(4, 1, 4)
        %     plot(time_trig + time_1st_trig, areas_pupil(:, 1)*1e1);
        %     hold on
        %     plot(time_trig + time_1st_trig, areas_pupil(:, 1)*1e1,'.', 'MarkerSize',15);
        %     title('DLC pupil area_x')
        %
        %     ax = findobj(gcf, 'type', 'axes');
        %     linkaxes(ax, 'x');
        % xlim([0 100])
        % xlim([5 12])
        %     xlim([10912 10916]) %pour 0208
        %     xlim([10630 10645]) %pour 0303
        %     xlim([5395 5400]) %pour 0508_2
        
    end
    
    % disp([num2str(time_1st_trig/1e4) ' s; nframes d: ' num2str(length(time_trig)-length(data_csv))])
    
    
end
% Save LFP-derived timing data for later use.
if exist(fullfile(datapath, 'video', 'DLC_data.mat'), 'file')
    save(fullfile(datapath, 'video', 'DLC_data.mat'), 'time_1st_trig', 'time_trig', '-append');
else
    save(fullfile(datapath, 'video', 'DLC_data.mat'), 'time_1st_trig', 'time_trig');
end

end