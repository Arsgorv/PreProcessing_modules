function Master_LFP_NP_preproc(sessions, opts)
% Master_LFP_NP_preproc
% Export Neuropixels LFP (Probe*-LFP) onto the MASTER (Acquisition Board) stitched timeline.
%
% Outputs:
%   per_channel:
%     <datapath>/ephys/LFPDataNP/LFP<ch>.mat      (variable LFP as tsd)
%   matrix (recommended for many channels):
%     <datapath>/ephys/LFPDataNP/LFPmat_<seg>.mat (variables Y int16, t0_ts, dt_ts, nOut, chUse, Fs_ds)
%
% Required:
%   opts.np_channels : 1-based channels within Probe*-LFP stream
%
% Optional:
%   opts.lfp_fs      : target Fs (default 1250)
%   opts.force       : overwrite existing (default false)
%   opts.maxChunkMB  : fread buffer budget (default 128)
%   opts.chunk_s     : cap chunk duration seconds (default 60)
%   opts.save_style  : 'per_channel' | 'matrix' | 'auto' (default 'per_channel')
%   opts.auto_threshold : if 'auto', >= threshold => matrix (default 64)
%   opts.allow_large_per_channel : allow per_channel for many chans (default false)

if nargin < 2, opts = struct(); end
if ischar(sessions) || isstring(sessions), sessions = {char(sessions)}; end

if ~isfield(opts,'np_channels') || isempty(opts.np_channels)
    error('Master_LFP_NP_preproc:NeedChannels','opts.np_channels is required.');
end
if ~isfield(opts,'lfp_fs'),  opts.lfp_fs = 1250; end
if ~isfield(opts,'force'),   opts.force  = false; end
if ~isfield(opts,'maxChunkMB'), opts.maxChunkMB = 128; end
if ~isfield(opts,'chunk_s'), opts.chunk_s = 60; end
if ~isfield(opts,'save_style'), opts.save_style = 'per_channel'; end
if ~isfield(opts,'auto_threshold'), opts.auto_threshold = 64; end
if ~isfield(opts,'allow_large_per_channel'), opts.allow_large_per_channel = false; end
if ~isfield(opts,'probe'), opts.probe = 'auto'; end   % 'auto' | 'A' | 'B'
if ~isfield(opts,'segments'),      opts.segments      = {}; end      % optional: basenames of ephys subfolders, in concatenation order
if ~isfield(opts,'require_manifest'), opts.require_manifest = false; end % if true, error when RAExp manifest is missing
if ~isfield(opts,'align_to_master'), opts.align_to_master = true; end   % linear NP->master time warp per segment
if ~isfield(opts,'sync_method'),     opts.sync_method     = 'auto'; end  % 'auto' | 'ttl' | 'linear' | 'none'
if ~isfield(opts,'master_stream'),   opts.master_stream   = 'Acquisition_Board'; end
if ~isfield(opts,'sync_ttl_line'),   opts.sync_ttl_line   = []; end      % [] = any rising edge

TsRate = 1e4;

for s = 1:numel(sessions)
    datapath = sessions{s};
    disp('------------------------------------------')
    disp(['[Master_LFP_NP_preproc] ' datapath])
    
    outDir   = fullfile(datapath,'ephys','LFPDataNP');
    obDir    = fullfile(datapath,'ephys','LFPData');
    if ~exist(outDir,'dir'), mkdir(outDir); end
    % obDir is only written into when channel numbers don't collide with
    % existing OB LFPs. mkdir lazily inside save_lfp_dual.
    
    % ---------------- segments on MASTER axis ----------------
    segs = np_build_segments(datapath, opts);
    if isempty(segs)
        warning('[Master_LFP_NP_preproc] No segments to process for %s; skipping.', datapath);
        continue
    end
    fprintf('[Master_LFP_NP_preproc] %d segment(s) on master timeline:\n', numel(segs));
    for i = 1:numel(segs)
        fprintf('  %2d %-60s t0=%8.2fs  dur=%8.2fs\n', i, segs(i).name, segs(i).t0_s, segs(i).dur_s);
    end
    
    obLfp = dir(fullfile(datapath,'ephys','*.lfp'));
    if ~isempty(obLfp)
        obBytes = obLfp(1).bytes;
        % You'd need nChan_OB and Fs_OB from the OB .xml; skip if unavailable
        % obDur = obBytes / (2 * nChan_OB * Fs_OB);
        % mismatch = obDur - segs(end).t1_s;
        % warn if abs(mismatch) > 1
    end
    
    % ---------------- channel list ----------------
    chList = unique(opts.np_channels(:))';
    nSel = numel(chList);
    
    % ---------------- choose mode ----------------
    save_style = lower(opts.save_style);
    if strcmp(save_style,'auto')
        if nSel >= opts.auto_threshold
            save_style_use = 'matrix';
        else
            save_style_use = 'per_channel';
        end
    else
        save_style_use = save_style;
    end
    
    if strcmp(save_style_use,'per_channel') && (nSel >= opts.auto_threshold) && ~opts.allow_large_per_channel
        warning('Requested %d channels with per_channel. Switching to matrix mode. Set opts.allow_large_per_channel=true to override.', nSel);
        save_style_use = 'matrix';
    end
    
    Fs_ds_bySeg = nan(numel(segs),1);
    syncBySeg = nan(numel(segs),1);
    probeBySeg = cell(numel(segs),1);
    syncInfoBySeg = cell(numel(segs),1);
    matFiles = {};
    
    % =====================================================================
    % MATRIX MODE (stream-write to avoid huge RAM; safe for 384 channels)
    % =====================================================================
    if strcmp(save_style_use,'matrix')
        
        for i = 1:numel(segs)
            segFolder = fullfile(datapath,'ephys', segs(i).name);
            [streamRoot, probeUsed] = resolve_np_lfp_stream(segFolder, opts.probe);
            if isempty(streamRoot)
                warning('  no Probe%s-LFP stream in %s; skipping segment', opts.probe, segs(i).name);
                continue
            end
            
            [Fs, nCh] = oe_read_stream_info(streamRoot, NaN, NaN);
            if ~isfinite(Fs),  Fs  = 2500; end
            if ~isfinite(nCh), nCh = 384;  end
            
            % drift sanity print
            try
                npSamples   = oe_nSamples_stream(streamRoot, nCh);
                npDur_s     = double(npSamples) / double(Fs);
                masterDur_s = (segs(i).t1_ts - segs(i).t0_ts) / TsRate;
                drift_s     = npDur_s - masterDur_s;
                syncBySeg(i) = drift_s;
                if abs(drift_s) > 0.500
                    warning('[%s] NP vs master drift = %.3f s', segs(i).name, drift_s);
                end
            catch
                syncBySeg(i) = NaN;
            end
            
            R = max(1, round(Fs/opts.lfp_fs));
            Fs_ds = Fs/R;
            dt_ts = round(TsRate/Fs_ds);
            
            keep = (chList >= 1) & (chList <= nCh);
            if ~any(keep)
                warning('  none of requested channels exist in %s', segs(i).name);
                continue
            end
            chUse = chList(keep);
            
            outMat = fullfile(outDir, ['LFPmat_' segs(i).name '.mat']);
            if exist(outMat,'file')==2 && ~opts.force
                disp(['  skip existing mat: ' outMat]);
                matFiles{end+1} = outMat; continue %#ok<AGROW>
            end
            
            segStart_ts = double(segs(i).t0_ts);
            segStop_ts  = double(segs(i).t1_ts);
            
            % build warp once per segment
            wopts = struct('use_ttl_sync', ~strcmp(opts.sync_method,'linear') ...
                && ~strcmp(opts.sync_method,'none'), ...
                'sync_ttl_line', opts.sync_ttl_line);
            [warpFn, warpInfo] = np_master_clock_warp( ...
                segFolder, streamRoot, opts.master_stream, wopts);
            if strcmp(opts.sync_method,'none'), warpFn = @(t) t; warpInfo.method='identity'; end
            save_segment_sync(datapath, segs(i).name, warpFn, warpInfo);
            
            % channel-batched decimate + warp + write
            chBatch = max(1, min(numel(chUse), ...
                floor(opts.maxChunkMB * 1024^2 / (8 * (Fs * masterDur_s)))));
            % heuristic: keep peak RAM per batch ~ maxChunkMB worth of doubles
            [Fs_ds, dt_ts, nOut] = oe_write_lfp_matrix_warped( ...
                streamRoot, outMat, chUse, nCh, Fs, opts.lfp_fs, ...
                segStart_ts, segStop_ts, TsRate, warpFn, chBatch, ...
                opts.maxChunkMB, opts.chunk_s);
            
            t0_ts = segStart_ts;
            save(outMat, 't0_ts','dt_ts','nOut','chUse','Fs_ds','-append');
            
            Fs_ds_bySeg(i)   = Fs_ds;
            probeBySeg{i}    = probeUsed;
            syncInfoBySeg{i} = warpInfo;
            matFiles{end+1}  = outMat; %#ok<AGROW>
            disp(['  saved mat: ' outMat ' (sync method: ' warpInfo.method ')']);
        end
        
        % =====================================================================
        % PER-CHANNEL MODE (fast for ~10 channels; ensures no overlap)
        % =====================================================================
    elseif strcmp(save_style_use,'per_channel')
        
        parts = cell(1,nSel);
        for j = 1:nSel
            parts{j} = {};
        end
        
        for i = 1:numel(segs)
            segFolder = fullfile(datapath,'ephys', segs(i).name);
            [streamRoot, probeUsed] = resolve_np_lfp_stream(segFolder, opts.probe);
            if isempty(streamRoot)
                warning('  no Probe%s-LFP stream in %s; skipping segment', opts.probe, segs(i).name);
                continue
            end
            
            [Fs, nCh] = oe_read_stream_info(streamRoot, NaN, NaN);
            if ~isfinite(Fs),  Fs  = 2500; end
            if ~isfinite(nCh), nCh = 384;  end
            
            % ---- sync metric (NP stream vs master segment length) ----
            try
                npSamples   = oe_nSamples_stream(streamRoot, nCh);
                npDur_s     = double(npSamples) / double(Fs);
                masterDur_s = (segs(i).t1_ts - segs(i).t0_ts) / TsRate;
                drift_s     = npDur_s - masterDur_s;
                syncBySeg(i) = drift_s;
                if abs(drift_s) > 0.500
                    warning('[%s] NP LFP vs master drift = %.3f s (NP=%.3f, master=%.3f)', ...
                        segs(i).name, drift_s, npDur_s, masterDur_s);
                else
                    fprintf('  [%s] drift %+0.3f s (within tolerance)\n', segs(i).name, drift_s);
                end
            catch ME
                warning('[%s] could not measure drift: %s', segs(i).name, ME.message);
                syncBySeg(i) = NaN;
            end
            
            R     = max(1, round(Fs / opts.lfp_fs));
            Fs_ds = Fs / R;
            
            dt_ts_real = TsRate / Fs_ds;
            dt_ts      = round(dt_ts_real);
            if abs(dt_ts_real - dt_ts) > 1e-9
                warning('dt_ts not integer (TsRate/Fs_ds). TsRate=%g Fs_ds=%g dt_ts=%g', ...
                    TsRate, Fs_ds, dt_ts_real);
            end
            
            keep = (chList >= 1) & (chList <= nCh);
            if ~any(keep), continue; end
            chUse = chList(keep);
            
            segStart_ts = double(segs(i).t0_ts);
            segStop_ts  = double(segs(i).t1_ts);
            maxNOut     = floor((segStop_ts - segStart_ts) / dt_ts);
            if maxNOut < 1
                warning('Bad segment boundaries: %s (stop<=start)', segs(i).name);
                continue
            end
            
            % ---- read full segment for the requested channels (int16) ----
            Xs = oe_read_channels_full_segment( ...
                streamRoot, chUse, nCh, Fs, opts.maxChunkMB, opts.chunk_s, maxNOut*R);
            if isempty(Xs), continue; end
            
            % ---- zero-phase decimation (or passthrough if R==1) ----
            % ---- decimation ----
            if R > 1
                Y_np = zeros(ceil(size(Xs,2)/R), numel(chUse), 'double');
                for jj = 1:numel(chUse)
                    Y_np(:,jj) = decimate(double(Xs(jj,:)), R);
                end
            else
                Y_np = double(Xs.');
            end
            
            % ---- master-clock warp (TTL -> linear -> identity) ----
            wopts = struct('use_ttl_sync', ~strcmp(opts.sync_method,'linear') ...
                && ~strcmp(opts.sync_method,'none'), ...
                'sync_ttl_line', opts.sync_ttl_line);
            [warpFn, warpInfo] = np_master_clock_warp( ...
                segFolder, streamRoot, opts.master_stream, wopts);
            if strcmp(opts.sync_method,'none'), warpFn = @(t) t; warpInfo.method='identity'; end
            save_segment_sync(datapath, segs(i).name, warpFn, warpInfo);
            
            nMaster = floor((segStop_ts - segStart_ts) / dt_ts);
            t_np_s     = (0:size(Y_np,1)-1)' / Fs_ds;
            t_np_in_m  = warpFn(t_np_s);
            t_target_s = (0:nMaster-1)' / Fs_ds;
            
            Yseg = zeros(nMaster, numel(chUse), 'int16');
            for jj = 1:numel(chUse)
                yi = interp1(t_np_in_m, Y_np(:,jj), t_target_s, 'linear', 0);
                Yseg(:,jj) = int16(yi);
            end
            
            Fs_ds_bySeg(i) = Fs_ds;
            probeBySeg{i}  = probeUsed;
            syncInfoBySeg{i} = warpInfo; %#ok<AGROW>
            
            t_ts = segStart_ts + (0:size(Yseg,1)-1)' * dt_ts;
            
            
            % append per channel, trimming any accidental overlap
            for jj = 1:numel(chUse)
                ch = chUse(jj);
                j  = find(chList == ch, 1, 'first');
                y  = Yseg(:,jj);
                
                if ~isempty(parts{j})
                    prev = Range(parts{j}{end});
                    prevEnd = prev(end);
                    if t_ts(1) <= prevEnd
                        kk = (t_ts > prevEnd);
                        t_ts2 = t_ts(kk); y2 = y(kk);
                    else
                        t_ts2 = t_ts; y2 = y;
                    end
                else
                    t_ts2 = t_ts; y2 = y;
                end
                
                if ~isempty(t_ts2)
                    parts{j}{end+1} = tsd(t_ts2, y2);
                end
            end
        end
        
        % save per channel
        for j = 1:nSel
            ch = chList(j);
            if isempty(parts{j})
                warning('  no data exported for ch %d', ch);
                continue
            end
            
            LFP = parts{j}{1};
            if numel(parts{j}) > 1
                LFP = cat(parts{j}{:});
            end
            temp = double(Data(LFP));
            LFP = tsd(Range(LFP), temp);
            
            [writtenNP, writtenOB] = save_lfp_dual(outDir, obDir, ch, LFP, opts.force);
            if writtenOB
                disp(['  saved: LFPData/LFP'   num2str(ch) '.mat (OB-compat)']);
            end
            if writtenNP
                disp(['  saved: LFPDataNP/LFP' num2str(ch) '.mat']);
            end
        end
        
    else
        error('Unknown opts.save_style: %s', opts.save_style);
    end
    
    % ---------------- Info ----------------
    Info = struct();
    Info.np_channels = chList;
    Info.lfp_fs = opts.lfp_fs;
    Info.maxChunkMB = opts.maxChunkMB;
    Info.chunk_s = opts.chunk_s;
    Info.save_style_requested = opts.save_style;
    Info.save_style_used = save_style_use;
    Info.auto_threshold = opts.auto_threshold;
    Info.created = datestr(now);
    Info.Fs_ds_bySeg = Fs_ds_bySeg;
    Info.segments = {segs.name};
    Info.matFiles = matFiles;
    Info.probeBySeg = probeBySeg;
    Info.probe_requested = opts.probe;
    Info.sync_drift_s_bySeg = syncBySeg;
    Info.align_to_master = opts.align_to_master;
    Info.sync_method_bySeg = cellfun(@(c) safefield(c,'method','unknown'), syncInfoBySeg, 'UniformOutput', false);
    Info.sync_corr_bySeg   = cellfun(@(c) safefield(c,'peak_corr',NaN), syncInfoBySeg);
    Info.sync_n_evts_bySeg = cellfun(@(c) safefield(c,'n_events',0),    syncInfoBySeg);
    
    save(fullfile(outDir,'InfoLFP_NP.mat'), 'Info');
    
end

end

% ============================ helpers ============================
function v = safefield(s,f,default)
if isstruct(s) && isfield(s,f), v = s.(f); else, v = default; end
end

function [streamRoot, probeUsed] = resolve_np_lfp_stream(segFolder, probePref)
% Try requested probe first; fall back to the other probe if not present.
probeUsed = '';
streamRoot = '';
switch lower(probePref)
    case 'a', order = {'ProbeA-LFP','ProbeB-LFP'};
    case 'b', order = {'ProbeB-LFP','ProbeA-LFP'};
    otherwise, order = {'ProbeA-LFP','ProbeB-LFP'};  % 'auto'
end
for k = 1:numel(order)
    s = oe_find_stream(segFolder, order{k});
    if ~isempty(s)
        streamRoot = s;
        probeUsed  = order{k};
        return
    end
end
end

function streamRoot = oe_find_stream(segFolder, key)
% Find a stream folder under segFolder whose name contains `key` (case-insensitive
% substring). Returns the first match that contains a continuous.dat file.
% Tolerates both "fixed" (recording*/continuous/...) and Record-Node layouts.
streamRoot = '';
contDirs = enumerate_continuous_dirs(segFolder);
keyL = lower(key);

% Prefer prefix match (e.g. 'Acquisition_Board' -> 'Acquisition_Board-100.acquisition_board')
for i = 1:numel(contDirs)
    [~, nm, ext] = fileparts(contDirs{i});
    nmFull = [nm ext];
    if startsWith(lower(nmFull), keyL) && exist(fullfile(contDirs{i},'continuous.dat'),'file') == 2
        streamRoot = contDirs{i};
        return
    end
end
% Fall back to substring match anywhere in the folder name
for i = 1:numel(contDirs)
    [~, nm, ext] = fileparts(contDirs{i});
    nmFull = [nm ext];
    if contains(lower(nmFull), keyL) && exist(fullfile(contDirs{i},'continuous.dat'),'file') == 2
        streamRoot = contDirs{i};
        return
    end
end
end

function dirs = enumerate_continuous_dirs(segFolder)
% Return the list of leaf directories that sit directly under any .../continuous/
% folder for the given segment, across both layouts.
dirs = {};
% (a) "fixed" layout: <seg>/recording*/continuous/<stream>/
recs = dir(fullfile(segFolder,'recording*'));
recs = recs([recs.isdir] & ~ismember({recs.name},{'.','..'}));
for r = 1:numel(recs)
    contPath = fullfile(recs(r).folder, recs(r).name, 'continuous');
    if isfolder(contPath)
        sub = dir(contPath);
        sub = sub([sub.isdir] & ~ismember({sub.name},{'.','..'}));
        for s = 1:numel(sub)
            dirs{end+1} = fullfile(sub(s).folder, sub(s).name); %#ok<AGROW>
        end
    end
end
% (b) Record-Node layout: <seg>/Record Node */experiment*/recording*/continuous/<stream>/
nodes = dir(fullfile(segFolder,'Record Node*'));
nodes = nodes([nodes.isdir]);
for n = 1:numel(nodes)
    exps = dir(fullfile(nodes(n).folder, nodes(n).name, 'experiment*'));
    exps = exps([exps.isdir]);
    for e = 1:numel(exps)
        recs2 = dir(fullfile(exps(e).folder, exps(e).name, 'recording*'));
        recs2 = recs2([recs2.isdir]);
        for r = 1:numel(recs2)
            contPath = fullfile(recs2(r).folder, recs2(r).name, 'continuous');
            if isfolder(contPath)
                sub = dir(contPath);
                sub = sub([sub.isdir] & ~ismember({sub.name},{'.','..'}));
                for s = 1:numel(sub)
                    dirs{end+1} = fullfile(sub(s).folder, sub(s).name); %#ok<AGROW>
                end
            end
        end
    end
end
end

function [Fs, nCh] = oe_read_stream_info(streamRoot, FsFallback, nChFallback)
% Robust to stream folders containing dots: "OneBox-102.Probe*-LFP"
Fs = FsFallback;
nCh = nChFallback;

recDir = fileparts(fileparts(streamRoot)); % .../recordingX
oebin = fullfile(recDir,'structure.oebin');
if exist(oebin,'file') ~= 2, return; end

try
    txt = fileread(oebin);
catch
    return
end

streamFolder = regexp(streamRoot, '[^\\/]+$', 'match', 'once');
if isempty(streamFolder), return; end

dotix = find(streamFolder=='.', 1, 'last');
if ~isempty(dotix)
    streamShort = streamFolder(dotix+1:end);
else
    streamShort = streamFolder;
end

% Anchor by folder_name first (always present), then fall back to stream_name.
key1 = ['"folder_name": "' streamFolder '/"'];
key2 = ['"folder_name": "' streamFolder '"'];
key3 = ['"folder_name":"' streamFolder '/"'];
key4 = ['"folder_name":"' streamFolder '"'];
ix = [];
for cand = {key1,key2,key3,key4}
    ix = strfind(txt, cand{1});
    if ~isempty(ix), break; end
end
if isempty(ix)
    % legacy fallback: try several stream_name variants
    variants = {streamShort, ...
        strrep(streamShort,'_',' '), ...
        strrep(streamShort,'-',' '), ...
        lower(streamShort), upper(streamShort)};
    for v = 1:numel(variants)
        for q = ["""", "'"]                                             %#ok<NBRAK2>
            ix = strfind(txt, ['"stream_name": ' char(q) variants{v} char(q)]);
            if ~isempty(ix), break; end
            ix = strfind(txt, ['"stream_name":' char(q) variants{v} char(q)]);
            if ~isempty(ix), break; end
        end
        if ~isempty(ix), break; end
    end
    if isempty(ix), return; end
end
ix = ix(1);

w0 = max(1, ix-6000);
w1 = min(numel(txt), ix+6000);
win = txt(w0:w1);
pos = ix - w0 + 1;

[s0, ~, tok] = regexp(win, 'sample_rate\"?\s*:\s*([0-9]+\.?[0-9]*)', 'start','end','tokens');
if ~isempty(s0)
    k = find(s0 < pos, 1, 'last');
    if isempty(k), k = 1; end
    Fs = str2double(tok{k}{1});
end

[s1, ~, tok] = regexp(win, 'num_channels\"?\s*:\s*([0-9]+)', 'start','end','tokens');
if ~isempty(s1)
    k = find(s1 > pos, 1, 'first');
    if isempty(k), k = find(s1 < pos, 1, 'last'); end
    if isempty(k), k = 1; end
    nCh = str2double(tok{k}{1});
end
end

function segs = np_build_segments(datapath, opts)
% Build a segments list on the master timeline, in this priority order:
%   1. opts.segments                                   (explicit override)
%   2. analysis/run_manifest_RAexp.mat                 (RAExp pipeline)
%   3. ephys/ExpeInfo.mat FolderForConcatenation_Ephys (matches OB ndm_concatenate)
%   4. Interactive listdlg on ephys/ subfolders        (generic fallback)
%
% Returns a struct array with fields:
%   name, path, dur_s, t0_s, t1_s, t0_ts, t1_ts
% Durations are read from the Acquisition Board stream so NP LFP lands on
% exactly the same timeline as the OB .lfp concatenation.

TsRate = 1e4;
ephysDir = fullfile(datapath,'ephys');
segs = [];

% --- (1) explicit override ---
if isfield(opts,'segments') && ~isempty(opts.segments)
    names = opts.segments;
    if ischar(names), names = {names}; end
    segs = build_segs_from_names(ephysDir, names, TsRate);
    fprintf('[np_build_segments] using opts.segments (%d)\n', numel(segs));
    return
end

% --- (2) RAExp manifest ---
mf = fullfile(datapath,'analysis','run_manifest_RAexp.mat');
if exist(mf,'file') == 2
    Srun = load(mf,'RUN'); RUN = Srun.RUN;
    segs = RUN.segments;
    for i = 1:numel(segs)
        segs(i).t0_ts = segs(i).t0_s * TsRate;
        segs(i).t1_ts = segs(i).t1_s * TsRate;
        if ~isfield(segs,'path') || isempty(segs(i).path)
            segs(i).path = fullfile(ephysDir, segs(i).name);
        end
    end
    fprintf('[np_build_segments] using RAExp manifest (%d segments)\n', numel(segs));
    return
end

if isfield(opts,'require_manifest') && opts.require_manifest
    error('Master_LFP_NP_preproc:NeedManifest', ...
        'opts.require_manifest=true but no run_manifest_RAexp.mat at %s', mf);
end

% --- (3) ExpeInfo concatenation list ---
eiFile = fullfile(ephysDir,'ExpeInfo.mat');
if exist(eiFile,'file') == 2
    try
        S = load(eiFile,'ExpeInfo');
        folders = S.ExpeInfo.PreProcessingInfo.FolderForConcatenation_Ephys;
        if ~iscell(folders), folders = {folders}; end
        names = cell(1,numel(folders));
        ok = true(1,numel(folders));
        for i = 1:numel(folders)
            [~, names{i}] = fileparts(folders{i});
            ok(i) = exist(fullfile(ephysDir, names{i}), 'dir') == 7;
        end
        if any(ok) && all(ok)
            segs = build_segs_from_names(ephysDir, names, TsRate);
            fprintf('[np_build_segments] using ExpeInfo.FolderForConcatenation_Ephys (%d segments)\n', numel(segs));
            return
        else
            fprintf('[np_build_segments] ExpeInfo list does not match ephys/ subfolders; falling back to picker\n');
            for i = 1:numel(folders)
                if ~ok(i), fprintf('    missing: %s\n', names{i}); end
            end
        end
    catch ME
        warning('Could not parse ExpeInfo for segment list: %s', ME.message);
    end
end

% --- (4) interactive picker ---
segs = pick_segments_interactive(ephysDir, TsRate);
end

function Xs = oe_read_channels_full_segment(streamRoot, chList, nCh, Fs, maxChunkMB, chunk_s, maxSamples)
% Read the full segment for a subset of channels into memory as int16.
% Xs is [nSel x nSamples]. Uses chunked fread so total peak RAM ~= chunk size.
Xs = [];
datFile = fullfile(streamRoot,'continuous.dat');
if exist(datFile,'file') ~= 2, return; end

info = dir(datFile);
nSamples = floor(double(info.bytes) / (2*double(nCh)));
if isfinite(maxSamples) && maxSamples > 0
    nSamples = min(nSamples, maxSamples);
end
if nSamples <= 0, return; end

chList = unique(chList(:))';
nSel   = numel(chList);
if nSel == 0, return; end

[chunkN, ~] = pick_chunk_size(Fs, nCh, maxChunkMB, chunk_s, 1);

Xs = zeros(nSel, nSamples, 'int16');

fid = fopen(datFile,'r');
if fid < 0, Xs = []; return; end
cleaner = onCleanup(@() fclose(fid));  %#ok<NASGU>

pos = 1;
s0  = 1;
while s0 <= nSamples
    nRead = min(chunkN, nSamples - s0 + 1);
    X = fread(fid, [nCh nRead], 'int16=>int16');
    if isempty(X), break; end
    k = size(X,2);
    Xs(:, pos:pos+k-1) = X(chList, :);
    pos = pos + k;
    s0  = s0 + nRead;
end
Xs = Xs(:, 1:pos-1);
end

function segs = build_segs_from_names(ephysDir, names, TsRate)
segs = struct('name',{},'path',{},'dur_s',{},'t0_s',{},'t1_s',{},'t0_ts',{},'t1_ts',{});
t0 = 0;
for i = 1:numel(names)
    nm = names{i};
    segPath = fullfile(ephysDir, nm);
    if exist(segPath,'dir') ~= 7
        warning('Segment folder not found, skipping: %s', segPath);
        continue
    end
    dur_s = oe_segment_acq_duration(segPath);
    if ~isfinite(dur_s) || dur_s <= 0
        warning('Could not determine duration for %s; skipping', nm);
        continue
    end
    segs(end+1).name  = nm;                                 %#ok<AGROW>
    segs(end).path    = segPath;
    segs(end).dur_s   = dur_s;
    segs(end).t0_s    = t0;
    segs(end).t1_s    = t0 + dur_s;
    segs(end).t0_ts   = segs(end).t0_s * TsRate;
    segs(end).t1_ts   = segs(end).t1_s * TsRate;
    t0 = segs(end).t1_s;
end
end

function segs = pick_segments_interactive(ephysDir, TsRate)
D = dir(ephysDir);
D = D([D.isdir]);
D = D(~ismember({D.name}, ...
    {'.','..','LFPData','LFPDataNP','ChannelsToAnalyse','analysis','figures','video'}));
if isempty(D)
    error('np_build_segments:NoSubfolders','No candidate subfolders under %s', ephysDir);
end

% Sort by parsed 'yyyy-mm-dd_HH-MM-SS' if present, else by name
tSort = nan(numel(D),1);
for i = 1:numel(D)
    tok = regexp(D(i).name, '(\d{4}-\d{2}-\d{2})_(\d{2}-\d{2}-\d{2})', 'tokens','once');
    if ~isempty(tok)
        try
            tSort(i) = datenum([tok{1} ' ' strrep(tok{2},'-',':')],'yyyy-mm-dd HH:MM:SS'); %#ok<DATNM>
        catch
        end
    end
end
if all(~isfinite(tSort))
    [~,ord] = sort({D.name});
else
    [~,ord] = sort(tSort);
end
D = D(ord);
names = {D.name};

[sel, ok] = listdlg( ...
    'PromptString', {'Select ephys segments to concatenate', ...
    'Hold Ctrl/Shift to pick multiple. Order = list order.'}, ...
    'ListString', names, ...
    'SelectionMode','multiple', ...
    'InitialValue', 1:numel(names), ...
    'ListSize',[560 420], ...
    'Name','NP LFP: select segments');
if ~ok || isempty(sel)
    error('np_build_segments:UserCancel','User cancelled segment selection.');
end
segs = build_segs_from_names(ephysDir, names(sel), TsRate);
end

function dur_s = oe_segment_acq_duration(segPath)
% Prefer Acquisition Board duration (so timeline matches OB .lfp concatenation),
% then fall back to probe streams if AcqBoard is missing.
keys = {'Acquisition_Board','Rhythm_FPGA','ProbeA-LFP','ProbeB-LFP'};
dur_s = NaN;
for k = 1:numel(keys)
    streamRoot = oe_find_stream(segPath, keys{k});
    if isempty(streamRoot), continue; end
    [Fs, nCh] = oe_read_stream_info(streamRoot, NaN, NaN);
    if ~isfinite(Fs), continue; end
    if ~isfinite(nCh), nCh = 80; end
    nSamp = oe_nSamples_stream(streamRoot, nCh);
    if isfinite(nSamp) && nSamp > 0
        dur_s = double(nSamp) / double(Fs);
        return
    end
end
end

function dur_s = oe_segment_duration_seconds(segFolder, streamKey)
streamRoot = oe_find_stream(segFolder, streamKey);
if isempty(streamRoot)
    dur_s = NaN; return
end

[Fs, nCh] = oe_read_stream_info(streamRoot, NaN, NaN);
if ~isfinite(Fs), Fs = 30000; end
if ~isfinite(nCh), nCh = 80; end

nSamp = oe_nSamples_stream(streamRoot, nCh);
dur_s = double(nSamp) / double(Fs);
end

function n = oe_nSamples_stream(streamRoot, nCh)
ts = dir(fullfile(streamRoot,'timestamps.npy'));
if ~isempty(ts)
    n = npy_numel(fullfile(ts(1).folder, ts(1).name));
    return
end
datFile = fullfile(streamRoot,'continuous.dat');
info = dir(datFile);
n = floor(double(info.bytes) / (2*double(nCh)));
end

function n = npy_numel(npyPath)
fid = fopen(npyPath,'r');
if fid < 0, error('Cannot open %s', npyPath); end
magic = fread(fid,6,'uint8=>char')';
if ~strcmp(magic, char([147 'NUMPY'])), fclose(fid); error('Bad npy: %s', npyPath); end
ver = fread(fid,2,'uint8');
if ver(1)==1
    hlen = fread(fid,1,'uint16');
else
    hlen = fread(fid,1,'uint32');
end
hdr = fread(fid, double(hlen), 'uint8=>char')';
fclose(fid);
tok = regexp(hdr, 'shape\W*\(\s*([0-9]+)', 'tokens','once');
n = str2double(tok{1});
end

function [Y_ds, Fs_ds] = oe_extract_channels_lfp_fread(streamRoot, chList, nCh, Fs, targetFs, maxChunkMB, chunk_s, maxNOut)
% Chunked extraction with anti-alias lowpass before decimation.
% If maxNOut provided, do not exceed master segment boundary.

datFile = fullfile(streamRoot,'continuous.dat');
if exist(datFile,'file')~=2
    Y_ds = []; Fs_ds = NaN; return
end

if nargin < 8, maxNOut = NaN; end
if nargin < 7 || isempty(chunk_s), chunk_s = 60; end
if nargin < 6 || isempty(maxChunkMB), maxChunkMB = 128; end
if nargin < 5 || isempty(targetFs), targetFs = Fs; end

chList = unique(chList(:))';
nSel = numel(chList);
if nSel == 0
    Y_ds = []; Fs_ds = NaN; return
end

R = max(1, round(Fs/targetFs));
Fs_ds = Fs / R;

info = dir(datFile);
nSamples = floor(double(info.bytes) / (2*double(nCh)));
if nSamples <= 0
    Y_ds = []; return
end

if isfinite(maxNOut) && maxNOut > 0
    nSamples = min(nSamples, maxNOut * R);
end

[chunkN, R] = pick_chunk_size(Fs, nCh, maxChunkMB, chunk_s, R);

% Anti-alias filter (active only when R > 1)
[bAA, grpDelay] = design_aa_fir(R);
haveAA = ~isempty(bAA);
if haveAA
    aaState = zeros(numel(bAA)-1, nSel);
else
    aaState = [];
end

% Over-allocate outputs by grpDelay/R samples; we will trim group delay at end.
nOutEst = floor(nSamples / R) + 8;
Y_ds = zeros(nOutEst, nSel, 'int16');

fid = fopen(datFile,'r');
if fid < 0, Y_ds = []; Fs_ds = NaN; return; end
cleaner = onCleanup(@() fclose(fid));

pos = 1;
s0 = 1;
phaseOffset = 0;   % sample index within input for next kept output
while s0 <= nSamples
    nRead = min(chunkN, nSamples - s0 + 1);
    X = fread(fid, [nCh nRead], 'int16=>int16');
    if isempty(X), break; end
    
    Xs = X(chList, :);                       % [nSel x nRead]
    if haveAA
        [Xf, aaState] = filter(bAA, 1, double(Xs), aaState, 2);  % filter along time
    else
        Xf = double(Xs);
    end
    
    % Sub-sample with continuous phase across chunk boundaries
    idx = (1 + phaseOffset):R:size(Xf,2);
    if ~isempty(idx)
        Xd = int16(Xf(:, idx));
        k = size(Xd, 2);
        Y_ds(pos:pos+k-1, :) = Xd.';
        pos = pos + k;
    end
    
    % Update phase for next chunk
    usedUp = size(Xf,2) - (isempty(idx)*0 + (~isempty(idx) * idx(end)));
    phaseOffset = mod(R - usedUp - 1, R);     % ensure next kept sample is R after last
    
    s0 = s0 + nRead;
end

Y_ds = Y_ds(1:pos-1, :);

% Compensate filter group delay: drop grpDelay/R samples from the start
if haveAA
    shift = round(grpDelay / R);
    if shift > 0 && size(Y_ds,1) > shift
        Y_ds = Y_ds(shift+1:end, :);
    end
end

% Final clamp to master segment budget
if isfinite(maxNOut) && maxNOut > 0 && size(Y_ds,1) > maxNOut
    Y_ds = Y_ds(1:maxNOut, :);
end
end

function [Fs_ds, dt_ts, nOut] = oe_write_lfp_matrix_fread(streamRoot, outMat, chList, nCh, Fs, targetFs, maxChunkMB, chunk_s, segStart_ts, segStop_ts, TsRate)
% Stream-write filtered+decimated Y into outMat, clamped to master segment stop.

datFile = fullfile(streamRoot,'continuous.dat');
info = dir(datFile);
nSamplesFile = floor(double(info.bytes) / (2*double(nCh)));

R = max(1, round(Fs/targetFs));
Fs_ds = Fs / R;

dt_ts_real = TsRate / Fs_ds;
dt_ts = round(dt_ts_real);
if abs(dt_ts_real - dt_ts) > 1e-9
    warning('dt_ts not integer (TsRate/Fs_ds). TsRate=%g Fs_ds=%g dt_ts=%g', TsRate, Fs_ds, dt_ts_real);
end

maxNOut = floor((segStop_ts - segStart_ts) / dt_ts);
if maxNOut < 1
    error('oe_write_lfp_matrix_fread:BadBounds','stop<=start');
end

nSamplesLimit = min(nSamplesFile, maxNOut * R);
nOutEst = floor(nSamplesLimit / R) + 8;
nSel = numel(chList);

M = matfile(outMat, 'Writable', true);
M.Y = zeros(nOutEst, nSel, 'int16');

[chunkN, R] = pick_chunk_size(Fs, nCh, maxChunkMB, chunk_s, R);

[bAA, grpDelay] = design_aa_fir(R);
haveAA = ~isempty(bAA);
if haveAA
    aaState = zeros(numel(bAA)-1, nSel);
else
    aaState = [];
end

fid = fopen(datFile,'r');
if fid < 0, error('oe_write_lfp_matrix_fread:OpenFail','Cannot open %s', datFile); end
cleaner = onCleanup(@() fclose(fid));

pos = 1;
s0 = 1;
phaseOffset = 0;
while s0 <= nSamplesLimit
    nRead = min(chunkN, nSamplesLimit - s0 + 1);
    X = fread(fid, [nCh nRead], 'int16=>int16');
    if isempty(X), break; end
    
    Xs = X(chList, :);
    if haveAA
        [Xf, aaState] = filter(bAA, 1, double(Xs), aaState, 2);
    else
        Xf = double(Xs);
    end
    
    idx = (1 + phaseOffset):R:size(Xf,2);
    if ~isempty(idx)
        Xd = int16(Xf(:, idx));
        k = size(Xd,2);
        M.Y(pos:pos+k-1, :) = Xd.';
        pos = pos + k;
    end
    
    if ~isempty(idx)
        phaseOffset = mod(R - (size(Xf,2) - idx(end)) - 1, R);
    end
    s0 = s0 + nRead;
end

% Group-delay is NOT trimmed: the linear-phase FIR introduces a constant
% delay of grpDelay/Fs seconds in the output time axis. Downstream code
% should subtract this offset from t0_ts if ms-level alignment matters.
M.filter_group_delay_s = grpDelay / Fs;

nOut = pos - 1;

if nOut > maxNOut
    Y = M.Y;
    Y = Y(1:maxNOut, :);
    nOut = maxNOut;
    M.Y = Y;
end
end

function [Fs_ds, dt_ts, nOut] = oe_write_lfp_matrix_warped( ...
    streamRoot, outMat, chList, nCh, Fs, targetFs, segStart_ts, segStop_ts, ...
    TsRate, warpFn, chBatchSize, maxChunkMB, chunk_s)
% Channel-batched matrix mode with master-clock warp.
% Reads NP raw, decimates per channel, applies warpFn, writes int16 column.

R = max(1, round(Fs/targetFs));
Fs_ds = Fs/R;
dt_ts_real = TsRate/Fs_ds;
dt_ts = round(dt_ts_real);

masterDur_s = (segStop_ts - segStart_ts) / TsRate;
nMaster = floor(masterDur_s * Fs_ds);

nSel = numel(chList);
M = matfile(outMat, 'Writable', true);
M.Y = zeros(nMaster, nSel, 'int16');

t_target_s = (0:nMaster-1)' / Fs_ds;

for b0 = 1:chBatchSize:nSel
    b1 = min(b0 + chBatchSize - 1, nSel);
    chBatch = chList(b0:b1);

    Xs = oe_read_channels_full_segment(streamRoot, chBatch, nCh, Fs, ...
                                       maxChunkMB, chunk_s, NaN);
    if isempty(Xs), continue; end

    if R > 1
        Y_np = zeros(ceil(size(Xs,2)/R), numel(chBatch), 'double');
        for jj = 1:numel(chBatch)
            Y_np(:,jj) = decimate(double(Xs(jj,:)), R);
        end
    else
        Y_np = double(Xs.');
    end
    clear Xs

    t_np_s = (0:size(Y_np,1)-1)' / Fs_ds;
    t_np_in_m = warpFn(t_np_s);

    for jj = 1:numel(chBatch)
        yi = interp1(t_np_in_m, Y_np(:,jj), t_target_s, 'linear', 0);
        M.Y(:, b0 + jj - 1) = int16(yi);
    end
    fprintf('    batch %d-%d / %d done\n', b0, b1, nSel);
end

nOut = nMaster;
end

function save_segment_sync(datapath, segName, warpFn, warpInfo)
% Persist the warp anchors so spike times can be aligned later via np_apply_master_warp.
outDir = fullfile(datapath,'analysis');
if ~exist(outDir,'dir'), mkdir(outDir); end
outFile = fullfile(outDir, ['master_clock_sync_' segName '.mat']);

method = warpInfo.method;
S = struct();
S.method = method;
S.created = datestr(now);
switch method
    case {'ttl','oe-synchronizer'}
        S.np_anchor_s     = warpInfo.np_event_times_s(:);
        if isfield(warpInfo,'master_event_times_s')
            S.master_anchor_s = warpInfo.master_event_times_s(:);
        else
            % oe-synchronizer: regenerate dense anchors from warpFn for portability
            tt = linspace(0, warpInfo.npDur_s, 200)';
            S.np_anchor_s     = tt;
            S.master_anchor_s = warpFn(tt);
        end
        S.peak_corr = warpInfo.peak_corr;
        S.n_events  = warpInfo.n_events;
    case 'linear-ratio'
        S.warp_ratio  = warpInfo.masterDur_s / max(warpInfo.npDur_s, eps);
        S.npDur_s     = warpInfo.npDur_s;
        S.masterDur_s = warpInfo.masterDur_s;
    case {'identity-no-master','identity'}
        % nothing extra
    otherwise
        warning('save_segment_sync: unknown method %s', method);
end
save(outFile, '-struct', 'S');
end

function [bAA, grpDelay] = design_aa_fir(R)
% Return FIR anti-alias lowpass coefficients for integer decimation by R.
% Normalized cutoff at 0.45 * (1/R) of Nyquist (targetFs/2 * 0.9 absolute).
% Returns empty when R==1 (pass-through, no filter cost).
if R <= 1
    bAA = [];
    grpDelay = 0;
    return
end
N = max(30, 8*R);              % order scales with R; ~60 taps for R=2, 80 for R=10
N = N + mod(N,2);              % even -> odd number of taps -> linear phase, integer delay
Wn = 0.9 / R;                  % cutoff at 90% of new Nyquist
bAA = fir1(N, Wn);             % Hamming-windowed FIR
grpDelay = N/2;                % samples @ input Fs
end

function [chunkN, R] = pick_chunk_size(Fs, nCh, maxChunkMB, chunk_s, R)
maxBytes = maxChunkMB * 1024^2;
chunkN = floor(maxBytes / (2*double(nCh)));
chunkN = max(R, floor(chunkN/R)*R);
chunkN = max(chunkN, Fs * 2);
chunkN = min(chunkN, Fs * 120);
if nargin >= 4 && isfinite(chunk_s) && chunk_s > 0
    chunkN2 = floor(chunk_s * Fs);
    chunkN2 = max(R, floor(chunkN2/R)*R);
    chunkN = min(chunkN, chunkN2);
end
end

function [writtenNP, writtenOB] = save_lfp_dual(npDir, obDir, ch, LFP, force)
% Save NP LFP into LFPDataNP always. Also mirror into LFPData if no collision.
% Collision = LFPData/LFP<ch>.mat already exists on disk.
writtenNP = false;
writtenOB = false;

% NP copy (always)
if ~exist(npDir,'dir'), mkdir(npDir); end
npFile = fullfile(npDir, ['LFP' num2str(ch) '.mat']);
if exist(npFile,'file')~=2 || force
    save(npFile, 'LFP', '-v7.3');
    writtenNP = true;
end

% OB-compat copy, only if no collision
obFile = fullfile(obDir, ['LFP' num2str(ch) '.mat']);
if exist(obFile,'file') == 2
    % Collision: do not touch OB copy (even with force), leave OB world untouched.
    return
end
if ~exist(obDir,'dir'), mkdir(obDir); end
save(obFile, 'LFP', '-v7.3');
writtenOB = true;
end