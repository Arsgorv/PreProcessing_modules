function [warpFn, info] = np_master_clock_warp(segFolder, npStreamRoot, masterStreamKey, opts)
% np_master_clock_warp
% Build NP-clock -> master-clock time map for one segment, using the shared
% analog sync pulse (1 Hz TTL recorded on BOTH systems):
%   - MASTER : LFPData/LFP<cfg.OneBox>.mat   (Acquisition Board ADC pulse, already
%              extracted to LFP by the OB pipeline)
%   - OneBox : OneBox-ADC stream, channel opts.onebox_adc_idx (default 3 = "ADC2")
%
% This is the SAME algorithm RAE_fit_onebox_to_master_map.m uses for spikes, so
% NP LFP and NP spikes land on an identical master clock.
%
% Returns:
%   warpFn(t_np_s) -> t_master_s   : affine map  t_master = a + b*t_np  (seconds,
%                                    both relative to segment start)
%   info : struct with .method, .a_s, .b, .rmse_ms, .drift_ppm, .n_used,
%          .np_anchor_s, .master_anchor_s, .npDur_s, .masterDur_s
%
% Strategy (in order):
%   1. 'affine-onebox' : analog pulse fit (preferred, sub-ms).
%   2. 'linear-ratio'  : scale by measured duration ratio (no pulse available).
%   3. 'identity'      : last resort.
%
% IMPORTANT: this file is intentionally SELF-CONTAINED (its own OE I/O helpers)
% so it does not depend on subfunctions living inside Master_LFP_NP_preproc.m.
%
% Required opts (passed by Master_LFP_NP_preproc):
%   opts.segStart_ts, opts.segStop_ts : segment bounds on master timeline (1e4 ticks)
% Optional opts:
%   opts.sync_method     : 'auto'(default) | 'affine' | 'linear' | 'none'
%   opts.master_ttl_chan : LFP channel carrying the master pulse (overrides cfg.OneBox)
%   opts.onebox_adc_idx  : 1-based channel in OneBox-ADC stream (default 3)
%   opts.minISI_s        : min inter-pulse interval for edge dedup (default 0.4)
%   opts.shift_search    : DEPRECATED/unused. Pulses are now matched 1:1 by
%                          nearest-neighbour (anchored at simultaneous OE start)
%                          and warped piecewise; no integer-shift search is done.
%   opts.thr_master      : master edge threshold ([]=auto)
%   opts.thr_onebox      : OneBox edge threshold ([]=auto)
%   opts.chunk_s         : OneBox-ADC scan chunk seconds (default 20)
%   opts.maxChunkMB      : OneBox-ADC scan chunk budget MB (default 64)

if nargin < 4, opts = struct(); end
if ~isfield(opts,'sync_method'),    opts.sync_method    = 'auto'; end
if ~isfield(opts,'master_ttl_chan'),opts.master_ttl_chan= []; end
if ~isfield(opts,'onebox_adc_idx'), opts.onebox_adc_idx = 3; end
if ~isfield(opts,'minISI_s'),       opts.minISI_s       = 0.4; end
if ~isfield(opts,'shift_search'),   opts.shift_search   = -60:60; end  % widened: -5:5 can clip the true device-start offset
if ~isfield(opts,'thr_master'),     opts.thr_master     = []; end
if ~isfield(opts,'thr_onebox'),     opts.thr_onebox     = []; end
if ~isfield(opts,'chunk_s'),        opts.chunk_s        = 20; end
if ~isfield(opts,'maxChunkMB'),     opts.maxChunkMB     = 64; end
if ~isfield(opts,'segStart_ts'),    opts.segStart_ts    = []; end
if ~isfield(opts,'segStop_ts'),     opts.segStop_ts     = []; end

TsRate = 1e4;

% derive datapath / segName from segFolder = <datapath>/ephys/<segName>
[ephysDir, segName] = fileparts(segFolder);
datapath = fileparts(ephysDir);

info = struct('method','none','a_s',0,'b',1,'rmse_ms',NaN,'drift_ppm',NaN, ...
              'n_used',0,'np_anchor_s',[],'master_anchor_s',[], ...
              'npDur_s',NaN,'masterDur_s',NaN);

% ---- always compute durations (needed for fallback + reporting) ----
[Fs_np, nCh_np] = oe_read_stream_info(npStreamRoot, NaN, NaN);
if ~isfinite(Fs_np), Fs_np = 2500; end
nNP = oe_nSamples_stream(npStreamRoot, nCh_np);
info.npDur_s = double(nNP) / Fs_np;
if ~isempty(opts.segStart_ts) && ~isempty(opts.segStop_ts)
    info.masterDur_s = double(opts.segStop_ts - opts.segStart_ts) / TsRate;
else
    info.masterDur_s = info.npDur_s; % unknown bounds -> assume equal
end

% ---- explicit overrides ----
if strcmpi(opts.sync_method,'none')
    warpFn = @(t) t; info.method = 'identity'; return
end
if strcmpi(opts.sync_method,'linear')
    [warpFn, info] = linear_fallback(info); return
end

% ============================ 1) analog-pulse affine ============================
try
    % --- master TTL channel ---
    % Resolution order (no interactive prompts unless truly nothing is known):
    %   1. opts.master_ttl_chan  (explicit, preferred)
    %   2. saved <datapath>/channel_cfg.mat  (written by a previous get_trigger_config)
    %   3. get_trigger_config(datapath)  (may prompt)
    abChan = opts.master_ttl_chan;
    if isempty(abChan)
        cfgFile = fullfile(datapath,'channel_cfg.mat');
        if exist(cfgFile,'file') == 2
            try
                Scfg = load(cfgFile,'cfg');
                if isfield(Scfg,'cfg') && isfield(Scfg.cfg,'OneBox'), abChan = Scfg.cfg.OneBox; end
            catch
            end
        end
    end
    if isempty(abChan)
        cfg = get_trigger_config(datapath);
        abChan = cfg.OneBox;
    end
    if isstruct(abChan)
        if isfield(abChan,'chan'),       abChan = abChan.chan;
        elseif isfield(abChan,'channel'),abChan = abChan.channel;
        elseif isfield(abChan,'idx'),    abChan = abChan.idx; end
    end
    if isempty(abChan) || any(isnan(abChan))
        error('np_master_clock_warp:NoMasterChan','cfg.OneBox/master_ttl_chan not set');
    end
    abChan = double(abChan(1));

    ttlFile = fullfile(datapath,'ephys','LFPData', sprintf('LFP%d.mat', abChan));
    if exist(ttlFile,'file') ~= 2
        error('np_master_clock_warp:NoMasterTTL','Missing %s', ttlFile);
    end
    Sttl = load(ttlFile,'LFP'); Lttl = Sttl.LFP;

    % restrict to segment window if bounds were provided
    if ~isempty(opts.segStart_ts) && ~isempty(opts.segStop_ts)
        Lw = Restrict(Lttl, intervalSet(double(opts.segStart_ts), double(opts.segStop_ts)));
        seg0_ts = double(opts.segStart_ts);
    else
        Lw = Lttl;
        r0 = Range(Lttl); seg0_ts = r0(1);
    end
    xM     = double(Data(Lw));
    tM_ts  = double(Range(Lw));
    tM_rel = (tM_ts - seg0_ts) / TsRate;   % seconds, relative to segment start
    if isempty(xM), error('np_master_clock_warp:EmptyMaster','no master samples in window'); end

    thrM = opts.thr_master;
    if isempty(thrM), thrM = (prctile(xM,5) + prctile(xM,95))/2; end
    tEdgeM = detect_rising_edges_from_vector(tM_rel, xM, thrM, opts.minISI_s);

    % --- OneBox ADC pulses (segment-local stream) ---
    streamADC = oe_find_stream(segFolder, 'OneBox-ADC');
    if isempty(streamADC), error('np_master_clock_warp:NoOneBoxADC','no OneBox-ADC in %s', segFolder); end
    [FsA, nChA] = oe_read_stream_info(streamADC, NaN, NaN);
    if ~isfinite(FsA), FsA = 30000; end
    adcIdx = opts.onebox_adc_idx;
    if adcIdx < 1 || adcIdx > nChA
        error('np_master_clock_warp:BadADCIdx','ADC idx %d out of range 1..%d', adcIdx, nChA);
    end
    [tEdgeO, ~] = detect_onebox_adc_edges(streamADC, adcIdx, nChA, FsA, opts);

    % --- pair the SAME physical pulses on both clocks, then warp PIECEWISE ---
    % A regular 1 Hz pulse has no pattern to lock an integer shift onto, so
    % rmse-minimising shift search is degenerate. OpenEphys starts all streams
    % on one Record command, so pulse #k on master == pulse #k on OneBox.
    % We (1) index-pair the leading min(N) pulses for a robust coarse line,
    % (2) refine by nearest-neighbour to fix any dropped/extra edge, then
    % (3) warp by INTERPOLATING THROUGH the matched pulse pairs (piecewise
    %     linear) so slow nonlinear drift is absorbed, not just a global slope.
    tEdgeM = tEdgeM(:); tEdgeO = tEdgeO(:);
    if numel(tEdgeM) < 10 || numel(tEdgeO) < 10
        error('np_master_clock_warp:TooFewPulses','need >=10 pulses (M=%d, O=%d)', ...
            numel(tEdgeM), numel(tEdgeO));
    end

    % (1) coarse affine from index pairing (simultaneous-start assumption)
    K = min(numel(tEdgeM), numel(tEdgeO));
    p0 = polyfit(tEdgeO(1:K), tEdgeM(1:K), 1);   % tM ~ b0*tO + a0
    b0 = p0(1); a0 = p0(2);

    % (2) nearest-neighbour match every OneBox pulse to a master pulse
    predM = a0 + b0*tEdgeO;                       % predicted master time per OneBox pulse
    nnIdx = interp1(tEdgeM, (1:numel(tEdgeM))', predM, 'nearest', 'extrap');
    nnIdx = min(max(round(nnIdx),1), numel(tEdgeM));
    dt    = tEdgeM(nnIdx) - predM;
    tol_s = 0.25;                                 % must be < half the pulse period
    good  = abs(dt) < tol_s;
    oO = tEdgeO(good);                            % matched OneBox-clock pulse times
    oM = tEdgeM(nnIdx(good));                     % matched master-clock pulse times
    % keep strictly increasing (interp1 requires unique, sorted knots)
    [oO, iu] = unique(oO, 'stable'); oM = oM(iu);
    keep2 = [true; diff(oM) > 0]; oO = oO(keep2); oM = oM(keep2);
    if numel(oO) < 10, error('np_master_clock_warp:MatchFailed','too few matched pulses'); end

    % (3) piecewise-linear warp through the matched pulse pairs
    warpFn = @(t) interp1(oO, oM, t, 'linear', 'extrap');

    % diagnostics: global affine summary + residual of pairs around it
    pa = polyfit(oO, oM, 1); b = pa(1); a = pa(2);
    resAff = (a + b*oO) - oM;                     % residual vs single straight line
    info.method          = 'piecewise-onebox';
    info.a_s             = a;                      % affine summary (reference only)
    info.b               = b;
    info.drift_ppm       = (b - 1) * 1e6;
    info.affine_rmse_ms  = sqrt(mean(resAff.^2)) * 1000;   % nonlinearity left by a single line
    info.affine_maxabs_ms= max(abs(resAff)) * 1000;
    info.rmse_ms         = 0;                      % piecewise passes through every knot
    info.n_used          = numel(oO);
    info.master_ttl_chan = abChan;
    info.onebox_adc_idx  = adcIdx;
    info.np_anchor_s     = oO;                     % OneBox-rel pulse times (== NP-LFP-rel)
    info.master_anchor_s = oM;                     % master-rel pulse times
    fprintf(['  [sync %s] piecewise n=%d  affine b=%.9f drift=%.1f ppm  ' ...
             'single-line residual rmse=%.3f ms max=%.3f ms  a0=%.3fs\n'], ...
        segName, info.n_used, b, info.drift_ppm, info.affine_rmse_ms, info.affine_maxabs_ms, a);
    return
catch ME
    warning('np_master_clock_warp: analog sync failed (%s); falling back. (%s)', ME.identifier, ME.message);
end

% ============================ 2) linear-ratio fallback ============================
[warpFn, info] = linear_fallback(info);
end

% =============================================================================
function [warpFn, info] = linear_fallback(info)
if isfinite(info.npDur_s) && info.npDur_s > 0 && isfinite(info.masterDur_s)
    ratio = info.masterDur_s / info.npDur_s;
else
    ratio = 1;
end
warpFn = @(t) ratio * t;
info.method = 'linear-ratio';
info.a_s = 0; info.b = ratio;
info.drift_ppm = (ratio - 1) * 1e6;
end

% =============================================================================
function tEdge = detect_rising_edges_from_vector(t_s, x, thr, minISI_s)
above = x > thr;
ix = find(diff([0; above(:)]) == 1);
tEdge = t_s(ix);
if isempty(tEdge), return; end
keep = true(size(tEdge));
lastKept = 1;
for k = 2:numel(tEdge)
    if (tEdge(k) - tEdge(lastKept)) < minISI_s
        keep(k) = false;
    else
        lastKept = k;
    end
end
tEdge = tEdge(keep);
end

function [tEdgeO, thr] = detect_onebox_adc_edges(streamADC, adcIdx, nCh, Fs, opts)
datFile = fullfile(streamADC,'continuous.dat');
if exist(datFile,'file') ~= 2, error('Missing %s', datFile); end

maxBytes = opts.maxChunkMB * 1024^2;
chunkN = floor(maxBytes / (2*double(nCh)));
chunkN = max(1, chunkN);
chunkN = min(chunkN, floor(opts.chunk_s * Fs));

% threshold from first ~5 s
nInit = min(chunkN, floor(5*Fs));
fid = fopen(datFile,'r');
if fid < 0, error('Cannot open %s', datFile); end
X0 = fread(fid, [nCh nInit], 'int16=>double');
fclose(fid);
if isempty(X0), error('Empty OneBox ADC file'); end
x0 = X0(adcIdx,:);
thr = opts.thr_onebox;
if isempty(thr), thr = (prctile(x0,5) + prctile(x0,95))/2; end

info = dir(datFile);
nSamples = floor(double(info.bytes) / (2*double(nCh)));

fid = fopen(datFile,'r');
if fid < 0, error('Cannot open %s', datFile); end
prevAbove = false; lastEdgeSamp = -Inf;
minISI_samp = round(opts.minISI_s * Fs);
edges = zeros(1,1e5); nE = 0; s0 = 1;
while s0 <= nSamples
    nRead = min(chunkN, nSamples - s0 + 1);
    X = fread(fid, [nCh nRead], 'int16=>double');
    if isempty(X), break; end
    x = X(adcIdx,:);
    above = (x > thr);
    d = diff([prevAbove, above]);
    ix = find(d == 1);
    for k = 1:numel(ix)
        samp = (s0 - 1) + ix(k);
        if samp - lastEdgeSamp >= minISI_samp
            nE = nE + 1;
            if nE > numel(edges), edges = [edges, zeros(1,1e5)]; end %#ok<AGROW>
            edges(nE) = samp; lastEdgeSamp = samp;
        end
    end
    prevAbove = above(end);
    s0 = s0 + nRead;
end
fclose(fid);
edges = edges(1:nE);
tEdgeO = (double(edges) - 1) / double(Fs);   % seconds since segment start (OneBox clock)
end

function [tM2, tO2] = align_by_shift(tM, tO, sh)
if sh > 0
    tOa = tO(1+sh:end); tMa = tM;
elseif sh < 0
    tMa = tM(1-sh:end); tOa = tO;
else
    tMa = tM; tOa = tO;
end
n = min(numel(tMa), numel(tOa));
tM2 = tMa(1:n); tO2 = tOa(1:n);
end

% ===================== self-contained OE I/O helpers =====================
function streamRoot = oe_find_stream(segFolder, key)
streamRoot = '';
contDirs = enumerate_continuous_dirs(segFolder);
keyL = lower(key);
for i = 1:numel(contDirs)
    [~, nm, ext] = fileparts(contDirs{i});
    if startsWith(lower([nm ext]), keyL) && exist(fullfile(contDirs{i},'continuous.dat'),'file')==2
        streamRoot = contDirs{i}; return
    end
end
for i = 1:numel(contDirs)
    [~, nm, ext] = fileparts(contDirs{i});
    if contains(lower([nm ext]), keyL) && exist(fullfile(contDirs{i},'continuous.dat'),'file')==2
        streamRoot = contDirs{i}; return
    end
end
end

function dirs = enumerate_continuous_dirs(segFolder)
dirs = {};
recs = dir(fullfile(segFolder,'recording*'));
recs = recs([recs.isdir] & ~ismember({recs.name},{'.','..'}));
for r = 1:numel(recs)
    contPath = fullfile(recs(r).folder, recs(r).name, 'continuous');
    if isfolder(contPath)
        sub = dir(contPath); sub = sub([sub.isdir] & ~ismember({sub.name},{'.','..'}));
        for s = 1:numel(sub), dirs{end+1} = fullfile(sub(s).folder, sub(s).name); end %#ok<AGROW>
    end
end
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
                sub = dir(contPath); sub = sub([sub.isdir] & ~ismember({sub.name},{'.','..'}));
                for s = 1:numel(sub), dirs{end+1} = fullfile(sub(s).folder, sub(s).name); end %#ok<AGROW>
            end
        end
    end
end
end

function [Fs, nCh] = oe_read_stream_info(streamRoot, FsFallback, nChFallback)
Fs = FsFallback; nCh = nChFallback;
recDir = fileparts(fileparts(streamRoot));
oebin = fullfile(recDir,'structure.oebin');
if exist(oebin,'file') ~= 2, return; end
try, txt = fileread(oebin); catch, return; end

streamFolder = regexp(streamRoot, '[^\\/]+$', 'match', 'once');
if isempty(streamFolder), return; end
dotix = find(streamFolder=='.', 1, 'last');
if ~isempty(dotix), streamShort = streamFolder(dotix+1:end); else, streamShort = streamFolder; end

% Anchor on folder_name first (always present), then stream_name variants.
ix = [];
for cand = {['"folder_name": "' streamFolder '/"'], ['"folder_name": "' streamFolder '"'], ...
            ['"folder_name":"'  streamFolder '/"'], ['"folder_name":"'  streamFolder '"']}
    ix = strfind(txt, cand{1});
    if ~isempty(ix), break; end
end
if isempty(ix)
    variants = {streamShort, strrep(streamShort,'_',' '), strrep(streamShort,'-',' '), ...
                lower(streamShort), upper(streamShort)};
    for v = 1:numel(variants)
        ix = strfind(txt, ['"stream_name": "' variants{v} '"']);
        if isempty(ix), ix = strfind(txt, ['"stream_name":"' variants{v} '"']); end
        if ~isempty(ix), break; end
    end
    if isempty(ix), return; end
end
ix = ix(1);
w0 = max(1, ix-6000); w1 = min(numel(txt), ix+6000);
win = txt(w0:w1); pos = ix - w0 + 1;

[s0, ~, tok] = regexp(win, 'sample_rate\"?\s*:\s*([0-9]+\.?[0-9]*)', 'start','end','tokens');
if ~isempty(s0)
    k = find(s0 < pos, 1, 'last'); if isempty(k), k = 1; end
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

function n = oe_nSamples_stream(streamRoot, nCh)
ts = dir(fullfile(streamRoot,'timestamps.npy'));
if ~isempty(ts)
    try
        n = local_npy_numel(fullfile(ts(1).folder, ts(1).name)); return
    catch
    end
end
datFile = fullfile(streamRoot,'continuous.dat');
info = dir(datFile);
if isempty(info) || ~isfinite(nCh) || nCh < 1, n = 0; return; end
n = floor(double(info.bytes) / (2*double(nCh)));
end

function n = local_npy_numel(npyPath)
fid = fopen(npyPath,'r');
if fid < 0, error('Cannot open %s', npyPath); end
magic = fread(fid,6,'uint8=>char')';
if ~strcmp(magic, char([147 'NUMPY'])), fclose(fid); error('Bad npy: %s', npyPath); end
ver = fread(fid,2,'uint8');
if ver(1)==1, hlen = fread(fid,1,'uint16'); else, hlen = fread(fid,1,'uint32'); end
hdr = fread(fid, double(hlen), 'uint8=>char')';
fclose(fid);
tok = regexp(hdr, 'shape\W*\(\s*([0-9]+)', 'tokens','once');
n = str2double(tok{1});
end
