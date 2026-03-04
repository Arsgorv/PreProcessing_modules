function OUT = RAE_apply_map_to_kilosort(datapath, segName, ksDir, mapFile, opts)
% RAE_apply_map_to_kilosort
% Convert Kilosort spikes (OneBox AP sample indices) to MASTER timeline using MAP from
% RAE_fit_onebox_to_master_map.
%
% Inputs:
%   ksDir   : folder containing spike_times.npy, spike_clusters.npy (Kilosort)
%   mapFile : .mat with struct MAP (a_s, b, segStart_ts)
%
% Optional:
%   opts.fs_ap : AP sampling rate for spike_times (default 30000)
%
% Output:
%   OUT struct saved to <datapath>/ephys/NP_spikes/synced/NP_synced_<segName>.mat

if nargin < 5, opts = struct(); end
if ~isfield(opts,'fs_ap'), opts.fs_ap = 30000; end

if exist(mapFile,'file') ~= 2
    error('RAE_apply_map_to_kilosort:NoMap','Missing mapFile: %s', mapFile);
end
S = load(mapFile,'MAP');
MAP = S.MAP;

if ~isfield(MAP,'a_s') || ~isfield(MAP,'b') || ~isfield(MAP,'segStart_ts')
    error('RAE_apply_map_to_kilosort:BadMap','MAP missing required fields.');
end

if exist(ksDir,'dir') ~= 7
    error('RAE_apply_map_to_kilosort:NoKSDir','Missing ksDir: %s', ksDir);
end

if exist('readNPY','file') ~= 2
    error('RAE_apply_map_to_kilosort:NoReadNPY','readNPY not found on path. Add npy toolbox / your packages.');
end

spikeTimesFile = fullfile(ksDir,'spike_times.npy');
spikeCluFile   = fullfile(ksDir,'spike_clusters.npy');
if exist(spikeTimesFile,'file') ~= 2
    error('Missing %s', spikeTimesFile);
end
if exist(spikeCluFile,'file') ~= 2
    error('Missing %s', spikeCluFile);
end

spk_samp = double(readNPY(spikeTimesFile));      % samples @ fs_ap
spk_clu  = double(readNPY(spikeCluFile));        % cluster ids

if numel(spk_samp) ~= numel(spk_clu)
    error('spike_times and spike_clusters length mismatch.');
end

t_onebox_rel = spk_samp / double(opts.fs_ap);              % seconds (OneBox clock, relative to segment start)
t_master_rel = MAP.a_s + MAP.b * t_onebox_rel;             % seconds on MASTER clock, relative to segment start
spike_ts = double(MAP.segStart_ts) + round(t_master_rel * 1e4);
spike_s  = spike_ts / 1e4;

OUT = struct();
OUT.datapath = datapath;
OUT.segName = segName;
OUT.ksDir = ksDir;
OUT.mapFile = mapFile;
OUT.fs_ap = opts.fs_ap;

OUT.map_a_s = MAP.a_s;
OUT.map_b = MAP.b;
OUT.map_shift = MAP.shift;
OUT.map_rmse_ms = MAP.rmse_ms;
OUT.map_drift_ppm = MAP.drift_ppm;

OUT.spike_ts = spike_ts;
OUT.spike_s = spike_s;
OUT.spike_clusters = spk_clu;

u = unique(spk_clu);
OUT.cluster_ids = u(:);
OUT.cluster_nspikes = zeros(numel(u),1);
for k = 1:numel(u)
    OUT.cluster_nspikes(k) = sum(spk_clu == u(k));
end

% optional: attach KS labels if present
OUT.cluster_info = [];
ci1 = fullfile(ksDir,'cluster_info.tsv');
ci2 = fullfile(ksDir,'cluster_group.tsv');
try
    if exist(ci1,'file')==2
        OUT.cluster_info = readtable(ci1,'FileType','text','Delimiter','\t');
    elseif exist(ci2,'file')==2
        OUT.cluster_info = readtable(ci2,'FileType','text','Delimiter','\t');
    end
catch
end

outDir = fullfile(datapath,'ephys','NP_spikes','synced');
if ~exist(outDir,'dir'), mkdir(outDir); end
outFile = fullfile(outDir, ['NP_synced_' segName '.mat']);
save(outFile,'OUT','-v7.3');

fprintf('[RAE_apply_map_to_kilosort] %s  spikes=%d  clusters=%d  saved=%s\n', ...
    segName, numel(spk_samp), numel(u), outFile);

end