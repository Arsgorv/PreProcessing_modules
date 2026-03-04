function [spikes, metadata] = load_wave_clus(directory)
% load_wave_clus
% Loads wave_clus outputs from:
%   directory/wave_clus/times_*.mat
%   directory/wave_clus/<subfolder>/times_*.mat
%
% Returns:
%   spikes   : [maxNspk x Nunits] spike times (ms) with NaN padding
%   metadata : [Nunits x 2] [channel, cluster]

wcDir = fullfile(directory);
if ~exist(wcDir,'dir')
    error('load_wave_clus:NoWaveClus','Missing folder: %s', wcDir);
end

% gather candidate folders: wcDir itself + its subfolders (1 level)
folders = dir(wcDir);
folders = folders([folders.isdir]);
folders = folders(~ismember({folders.name},{'.','..'}));

searchFolders = cell(1, 1+numel(folders));
searchFolders{1} = wcDir;
for i = 1:numel(folders)
    searchFolders{1+i} = fullfile(wcDir, folders(i).name);
end

allSpikes = {};
metadata  = [];

for f = 1:numel(searchFolders)
    cur = searchFolders{f};
    timesFiles = dir(fullfile(cur, 'times_*.mat'));
    for t = 1:numel(timesFiles)
        fileName = timesFiles(t).name;
        match = regexp(fileName, 'times_C(\d+)\.mat', 'tokens','once');
        if isempty(match), continue; end
        ch = str2double(match{1});

        S = load(fullfile(timesFiles(t).folder, timesFiles(t).name), 'cluster_class');
        if ~isfield(S,'cluster_class'), continue; end
        cc = S.cluster_class;

        uClu = unique(cc(:,1));
        uClu(uClu == 0) = []; % noise cluster

        for c = 1:numel(uClu)
            clu = uClu(c);
            spk = cc(cc(:,1) == clu, 2); % ms (wave_clus convention)
            allSpikes{end+1} = spk(:); %#ok<AGROW>
            metadata = [metadata; ch, clu]; %#ok<AGROW>
        end
    end
end

if isempty(allSpikes)
    spikes = NaN(0,0);
    metadata = zeros(0,2);
    return
end

maxSpk = max(cellfun(@numel, allSpikes));
spikes = NaN(maxSpk, numel(allSpikes));
for i = 1:numel(allSpikes)
    spikes(1:numel(allSpikes{i}), i) = allSpikes{i};
end
end