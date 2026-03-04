function openEphys2wave_clus(streamRoot,chLst,numChannels,CommonRefChannels, outDir, Fs)
% openEphys2wave_clus
% Convert one OpenEphys "continuous stream" (folder with continuous.dat) into wave_clus input,
% run Get_spikes + Do_clustering per channel, and save outputs in outDir.
%
% streamRoot: .../recording1/continuous/Acquisition_Board-100.acquisition_board (folder containing continuous.dat)
% outDir    : folder where C*.mat + wave_clus outputs will be written

if nargin < 6 || isempty(Fs), Fs = 30000; end
if nargin < 5 || isempty(outDir), outDir = pwd; end
if nargin < 4, CommonRefChannels = []; end

if ~exist(outDir,'dir'), mkdir(outDir); end
cd(outDir);

datFile = fullfile(streamRoot, 'continuous.dat');
if ~exist(datFile,'file')
    error('openEphys2wave_clus:NoDat','Missing continuous.dat at %s', datFile);
end

bytesPerSample = 2; % int16
info = dir(datFile);
numSamples = info.bytes / (numChannels * bytesPerSample);
numSamples = floor(numSamples);

m = memmapfile(datFile, ...
    'Format', {'int16', [numChannels, numSamples], 'data'}, ...
    'Writable', false);

% --- compute common reference (optional) ---
segmentSize = Fs * 15; % 15s segments
if ~isempty(CommonRefChannels)
    disp('compute common ref...');
    commonRef = zeros(1, size(m.Data.data, 2), 'single');
    nSegments = ceil(size(m.Data.data, 2) / segmentSize);
    for iSeg = 1:nSegments
        s0 = (iSeg-1)*segmentSize + 1;
        s1 = min(iSeg*segmentSize, size(m.Data.data, 2));
        x = double(m.Data.data(CommonRefChannels, s0:s1));
        commonRef(s0:s1) = single(mean(x, 1));
    end
else
    commonRef = [];
end

% --- per channel ---
for cnum = chLst
    fprintf('\nSaving raw channel %d\n', cnum);

    x = double(m.Data.data(cnum, :));
    if ~isempty(commonRef)
        x = x - double(commonRef);
    end

    par = struct(); 
    sr  = Fs;
    data = x;

    rawdataSingleChannel = ['C' num2str(cnum) '.mat'];
    save(rawdataSingleChannel, 'data', 'sr', '-v7.3');

    fprintf('Processing channel %d\n', cnum);
    spikesSingleChannel = ['C' num2str(cnum) '_spikes.mat'];

    Get_spikes(rawdataSingleChannel, 'parallel', false);
    Do_clustering(spikesSingleChannel, 'parallel', false);
end

end