function [spikes, metadata] = load_wave_clus(directory)
    % Change to the wave_clus directory
    cd([directory '/wave_clus'])

    % Get a list of all folders in the current directory
    folders = dir;
    folders = folders([folders.isdir]); % Only keep directories
    folders = folders(~ismember({folders.name}, {'.', '..'})); % Exclude '.' and '..'

    % Initialize variables
    allSpikes = {}; % Temporary storage for spike times
    metadata = [];  % To map each column of spikes to channel and cluster

    % Loop through each folder
    for f = 1:length(folders)
        % Get the full path of the current folder
        currentFolder = fullfile(folders(f).folder, folders(f).name);

        % Find all files starting with "times_"
        timesFiles = dir(fullfile(currentFolder, 'times_*.mat'));

        % Loop through each time file
        for t = 1:length(timesFiles)
            % Extract the number from the filename (assuming "times_C*.mat")
            fileName = timesFiles(t).name;
            match = regexp(fileName, 'times_C(\d+)\.mat', 'tokens');
            if isempty(match)
                continue; % Skip if no match is found
            end
            ch = str2double(match{1}{1}); % Extract channel number as a numeric value

            % Load the cluster_class data
            clusterClassPath = fullfile(timesFiles(t).folder, timesFiles(t).name);
            disp(['loading ' timesFiles(t).name ' ...'])
            clusterData = load(clusterClassPath, 'cluster_class');
            clusterData = clusterData.cluster_class;

            % Find unique clusters
            uniqueClusters = unique(clusterData(:, 1));

            % Remove cluster #0
            uniqueClusters(uniqueClusters == 0) = [];

            % Store spikes data for each cluster
            for c = 1:length(uniqueClusters)
                clusterID = uniqueClusters(c);
                spikesInCluster = clusterData(clusterData(:, 1) == clusterID, 2); % Extract spike times

                % Append to temporary storage
                allSpikes{end+1} = spikesInCluster; % Append spike times
                metadata = [metadata; ch, clusterID]; % Store channel and cluster mapping
            end
        end
    end

    % Combine all spike times into a 2D matrix (padded with NaN)
    maxSpikes = max(cellfun(@length, allSpikes)); % Find the max number of spikes
    spikes = NaN(maxSpikes, length(allSpikes));  % Initialize with NaN
    for i = 1:length(allSpikes)
        spikes(1:length(allSpikes{i}), i) = allSpikes{i}; % Fill each column with spike times
    end
end
