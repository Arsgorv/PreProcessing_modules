function RES = Master_NP_spikes_sync_preproc(sessions, opts)
% Master_NP_spikes_sync_preproc
% For each session and each ephys segment (from run_manifest_RAexp.mat):
%   1) Fit OneBox->MASTER map using 1 Hz TTL (master: cfg.OneBox; onebox: ADC2)
%   2) Apply map to Kilosort spikes (spike_times.npy / spike_clusters.npy)
%
% Assumes Kilosort output is stored in:
%   <datapath>/ephys/NP_spikes/<segName>/
%
% Saves:
%   maps:    <datapath>/analysis/sync_onebox/map_<segName>.mat
%   spikes:  <datapath>/ephys/NP_spikes/synced/NP_synced_<segName>.mat
%   summary: <datapath>/ephys/NP_spikes/synced/NP_synced_summary.csv
%
% Required external:
%   - RAE_fit_onebox_to_master_map.m
%   - RAE_apply_map_to_kilosort.m
%   - readNPY.m (for .npy)

if nargin < 2, opts = struct(); end
if ischar(sessions) || isstring(sessions), sessions = {char(sessions)}; end

% ---- defaults ----
if ~isfield(opts,'onebox_adc_idx'), opts.onebox_adc_idx = 3; end   % ADC2 -> index 3
if ~isfield(opts,'fs_ap'),          opts.fs_ap = 30000; end
if ~isfield(opts,'force'),          opts.force = false; end
if ~isfield(opts,'shift_search'),   opts.shift_search = -5:5; end
if ~isfield(opts,'minISI_s'),       opts.minISI_s = 0.4; end
if ~isfield(opts,'chunk_s'),        opts.chunk_s = 20; end
if ~isfield(opts,'maxChunkMB'),     opts.maxChunkMB = 64; end
if ~isfield(opts,'ks_root'),        opts.ks_root = fullfile('ephys','NP_spikes'); end
if ~isfield(opts,'phase_filter'),   opts.phase_filter = {}; end  

if exist('RAE_fit_onebox_to_master_map','file')~=2
    error('Missing RAE_fit_onebox_to_master_map on path.');
end
if exist('RAE_apply_map_to_kilosort','file')~=2
    error('Missing RAE_apply_map_to_kilosort on path.');
end
if exist('readNPY','file')~=2
    error('Missing readNPY on path.');
end

rows = {};
bySeg = struct('datapath',{},'segName',{},'phase',{},'ksDir',{},'status',{}, ...
    'nSpikes',{},'nClusters',{},'rmse_ms',{},'maxabs_ms',{},'drift_ppm',{},'mapFile',{},'outFile',{});

for si = 1:numel(sessions)
    datapath = sessions{si};
    disp('------------------------------------------')
    disp(['[Master_NP_spikes_sync_preproc] ' datapath])

    mf = fullfile(datapath,'analysis','run_manifest_RAexp.mat');
    if exist(mf,'file')~=2
        warning('  missing run_manifest_RAexp.mat (skip): %s', datapath);
        continue
    end
    Srun = load(mf,'RUN'); RUN = Srun.RUN;
    segs = RUN.segments;

    outDir = fullfile(datapath,'ephys','NP_spikes','synced');
    if ~exist(outDir,'dir'), mkdir(outDir); end

    for i = 1:numel(segs)
        segName = segs(i).name;

        ph = '';
        if isfield(segs,'phase') && ~isempty(segs(i).phase)
            ph = segs(i).phase;
        end

        if ~isempty(opts.phase_filter)
            if isempty(ph) || ~any(strcmpi(ph, opts.phase_filter))
                continue
            end
        end

        ksDir = fullfile(datapath, opts.ks_root, segName);
        spkTimes = fullfile(ksDir,'spike_times.npy');
        spkClu   = fullfile(ksDir,'spike_clusters.npy');

        st = 'ok';
        mapFile = fullfile(datapath,'analysis','sync_onebox',['map_' segName '.mat']);
        outFile = fullfile(outDir, ['NP_synced_' segName '.mat']);

        if exist(ksDir,'dir')~=7 || exist(spkTimes,'file')~=2 || exist(spkClu,'file')~=2
            st = 'no_kilosort';
            rows(end+1,:) = {datapath, segName, ph, ksDir, st, NaN, NaN, NaN, NaN, NaN}; %#ok<AGROW>
            bySeg(end+1) = pack_bySeg(datapath, segName, ph, ksDir, st, mapFile, outFile); %#ok<AGROW>
            continue
        end

        % If outputs exist and not forcing, skip recompute but still summarize
        if exist(outFile,'file')==2 && exist(mapFile,'file')==2 && ~opts.force
            S = load(outFile,'OUT');
            OUT = S.OUT;
            rows(end+1,:) = {datapath, segName, ph, ksDir, 'exists', numel(OUT.spike_ts), numel(OUT.cluster_ids), ...
                OUT.map_rmse_ms, NaN, OUT.map_drift_ppm}; %#ok<AGROW>
            bySeg(end+1) = pack_bySeg(datapath, segName, ph, ksDir, 'exists', mapFile, outFile, OUT); %#ok<AGROW>
            continue
        end

        % ---- fit map ----
        mapOpts = struct();
        mapOpts.onebox_adc_idx = opts.onebox_adc_idx;
        mapOpts.shift_search   = opts.shift_search;
        mapOpts.minISI_s       = opts.minISI_s;
        mapOpts.chunk_s        = opts.chunk_s;
        mapOpts.maxChunkMB     = opts.maxChunkMB;

        try
            MAP = RAE_fit_onebox_to_master_map(datapath, segName, mapOpts); %#ok<NASGU>
        catch ME
            st = ['map_fail:' ME.identifier];
            rows(end+1,:) = {datapath, segName, ph, ksDir, st, NaN, NaN, NaN, NaN, NaN}; %#ok<AGROW>
            bySeg(end+1) = pack_bySeg(datapath, segName, ph, ksDir, st, mapFile, outFile); %#ok<AGROW>
            continue
        end

        % ---- apply to kilosort ----
        appOpts = struct();
        appOpts.fs_ap = opts.fs_ap;

        try
            OUT = RAE_apply_map_to_kilosort(datapath, segName, ksDir, mapFile, appOpts); %#ok<NASGU>
            % note: function saves OUT internally
            S = load(outFile,'OUT');
            OUT = S.OUT;

            rows(end+1,:) = {datapath, segName, ph, ksDir, 'ok', numel(OUT.spike_ts), numel(OUT.cluster_ids), ...
                OUT.map_rmse_ms, NaN, OUT.map_drift_ppm}; %#ok<AGROW>
            bySeg(end+1) = pack_bySeg(datapath, segName, ph, ksDir, 'ok', mapFile, outFile, OUT); %#ok<AGROW>
        catch ME
            st = ['apply_fail:' ME.identifier];
            rows(end+1,:) = {datapath, segName, ph, ksDir, st, NaN, NaN, NaN, NaN, NaN}; %#ok<AGROW>
            bySeg(end+1) = pack_bySeg(datapath, segName, ph, ksDir, st, mapFile, outFile); %#ok<AGROW>
            continue
        end
    end

    % per-session summary csv
    if ~isempty(rows)
        T = cell2table(rows, 'VariableNames', ...
            {'datapath','segment','phase','ksDir','status','nSpikes','nClusters','rmse_ms','maxabs_ms','drift_ppm'});
        writetable(T, fullfile(outDir,'NP_synced_summary.csv'));
    end
end

RES = struct();
if isempty(rows)
    RES.T = table();
else
    RES.T = cell2table(rows, 'VariableNames', ...
        {'datapath','segment','phase','ksDir','status','nSpikes','nClusters','rmse_ms','maxabs_ms','drift_ppm'});
end
RES.bySeg = bySeg;

end

% ---------------- helper ----------------
function S = pack_bySeg(datapath, segName, phase, ksDir, status, mapFile, outFile, OUT)
if nargin < 8, OUT = []; end
S = struct();
S.datapath = datapath;
S.segName = segName;
S.phase = phase;
S.ksDir = ksDir;
S.status = status;
S.mapFile = mapFile;
S.outFile = outFile;

S.nSpikes = NaN; S.nClusters = NaN; S.rmse_ms = NaN; S.maxabs_ms = NaN; S.drift_ppm = NaN;
if ~isempty(OUT)
    if isfield(OUT,'spike_ts'), S.nSpikes = numel(OUT.spike_ts); end
    if isfield(OUT,'cluster_ids'), S.nClusters = numel(OUT.cluster_ids); end
    if isfield(OUT,'map_rmse_ms'), S.rmse_ms = OUT.map_rmse_ms; end
    if isfield(OUT,'map_drift_ppm'), S.drift_ppm = OUT.map_drift_ppm; end
end
end