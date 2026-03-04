function MAP = RAE_fit_onebox_to_master_map(datapath, segName, opts)
% RAE_fit_onebox_to_master_map
% Fit OneBox time to MASTER time using 1 Hz TTL in:
%   - MASTER: Acquisition board LFP channel cfg.OneBox (stitched LFPData)
%   - OneBox: OneBox-ADC stream, ADC2 (default opts.onebox_adc_idx=3)
%
% Fits: t_master_rel = a + b * t_onebox_rel  (seconds, relative to seg start)
% Saves: <datapath>/analysis/sync_onebox/map_<segName>.mat

if nargin < 3, opts = struct(); end
if ~isfield(opts,'onebox_adc_idx'), opts.onebox_adc_idx = 3; end % ADC2 -> 3 (1-based)
if ~isfield(opts,'minISI_s'), opts.minISI_s = 0.4; end
if ~isfield(opts,'thr_master'), opts.thr_master = []; end
if ~isfield(opts,'thr_onebox'), opts.thr_onebox = []; end
if ~isfield(opts,'shift_search'), opts.shift_search = -5:5; end
if ~isfield(opts,'chunk_s'), opts.chunk_s = 20; end % for scanning OneBox ADC
if ~isfield(opts,'maxChunkMB'), opts.maxChunkMB = 64; end

TsRate = 1e4;

% --- load manifest + segment bounds ---
mf = fullfile(datapath,'analysis','run_manifest_RAexp.mat');
if exist(mf,'file') ~= 2
    error('RAE_fit_onebox_to_master_map:NeedManifest','Missing run_manifest_RAexp.mat');
end
Srun = load(mf,'RUN'); RUN = Srun.RUN;

segNames = {RUN.segments.name};
ix = find(strcmp(segNames, segName), 1, 'first');
if isempty(ix)
    error('RAE_fit_onebox_to_master_map:BadSeg','segName not found in manifest: %s', segName);
end
seg = RUN.segments(ix);

segStart_ts = seg.t0_s * TsRate;
segStop_ts  = seg.t1_s * TsRate;

% --- master TTL channel from config ---
cfg = get_trigger_config(datapath);
abChan = cfg.OneBox;
if isstruct(abChan)
    if isfield(abChan,'chan'), abChan = abChan.chan;
    elseif isfield(abChan,'channel'), abChan = abChan.channel;
    elseif isfield(abChan,'idx'), abChan = abChan.idx;
    else
        error('cfg.OneBox is a struct but no chan/channel/idx field found.');
    end
end
abChan = double(abChan);

ttlFile = fullfile(datapath,'ephys','LFPData', sprintf('LFP%d.mat', abChan));
if exist(ttlFile,'file') ~= 2
    error('RAE_fit_onebox_to_master_map:NoMasterTTL','Missing master TTL file: %s', ttlFile);
end
Sttl = load(ttlFile,'LFP');
Lttl = Sttl.LFP;

% restrict to segment
Lw = Restrict(Lttl, intervalSet(segStart_ts, segStop_ts));
xM = double(Data(Lw));
tM_ts = double(Range(Lw));
tM_rel = (tM_ts - segStart_ts) / TsRate; % seconds

if isempty(xM)
    error('RAE_fit_onebox_to_master_map:EmptyMaster','No master TTL samples in segment window.');
end

thrM = opts.thr_master;
if isempty(thrM)
    lo = prctile(xM, 5);
    hi = prctile(xM, 95);
    thrM = (lo + hi) / 2;
end

tEdgeM = detect_rising_edges_from_vector(tM_rel, xM, thrM, opts.minISI_s);

% --- OneBox ADC2 pulses (segment-local stream) ---
segFolder = fullfile(datapath,'ephys', segName);
streamADC = oe_find_stream(segFolder, 'OneBox-ADC');
if isempty(streamADC)
    error('RAE_fit_onebox_to_master_map:NoOneBoxADC','Cannot find OneBox-ADC stream in %s', segFolder);
end

[FsA, nChA] = oe_read_stream_info(streamADC, NaN, NaN);
if ~isfinite(FsA), FsA = 30000; end

adcIdx = opts.onebox_adc_idx;
if adcIdx < 1 || adcIdx > nChA
    error('ADC index %d out of range 1..%d', adcIdx, nChA);
end

[tEdgeO, thrO] = detect_onebox_adc_edges(streamADC, adcIdx, nChA, FsA, opts);

% --- align pulse trains by shift search and fit affine ---
best = struct('rmse',Inf,'a',NaN,'b',NaN,'shift',0,'nUsed',0,'tM',[],'tO',[]);
for sh = opts.shift_search
    [tM2, tO2] = align_by_shift(tEdgeM, tEdgeO, sh);
    if numel(tM2) < 10, continue; end

    p = polyfit(tO2(:), tM2(:), 1); % tM = b*tO + a
    b = p(1); a = p(2);
    pred = a + b*tO2(:);
    err = pred - tM2(:);
    rmse = sqrt(mean(err.^2));

    if rmse < best.rmse
        best.rmse = rmse;
        best.a = a;
        best.b = b;
        best.shift = sh;
        best.nUsed = numel(tM2);
        best.tM = tM2(:);
        best.tO = tO2(:);
        best.err = err;
    end
end

if ~isfinite(best.a)
    error('RAE_fit_onebox_to_master_map:FitFailed','Could not fit map (insufficient aligned pulses).');
end

MAP = struct();
MAP.datapath = datapath;
MAP.segName = segName;
MAP.segStart_ts = segStart_ts;
MAP.segStop_ts  = segStop_ts;
MAP.master_ttl_chan = abChan;
MAP.onebox_adc_idx = adcIdx;
MAP.Fs_onebox_adc = FsA;
MAP.thr_master = thrM;
MAP.thr_onebox = thrO;

MAP.a_s = best.a;
MAP.b = best.b;
MAP.shift = best.shift;
MAP.nUsed = best.nUsed;
MAP.rmse_ms = best.rmse * 1000;
MAP.maxabs_ms = max(abs(best.err)) * 1000;
MAP.drift_ppm = (MAP.b - 1) * 1e6;

MAP.tMaster_rel_used = best.tM;
MAP.tOneBox_rel_used = best.tO;
MAP.created = datestr(now);

outDir = fullfile(datapath,'analysis','sync_onebox');
if ~exist(outDir,'dir'), mkdir(outDir); end
outFile = fullfile(outDir, ['map_' segName '.mat']);
save(outFile,'MAP','-v7.3');

fprintf('[RAE_fit_onebox_to_master_map] %s  a=%.6f  b=%.9f  rmse=%.3f ms  drift=%.1f ppm  n=%d  shift=%d\n', ...
    segName, MAP.a_s, MAP.b, MAP.rmse_ms, MAP.drift_ppm, MAP.nUsed, MAP.shift);

end

% ================= helpers =================

function tEdge = detect_rising_edges_from_vector(t_s, x, thr, minISI_s)
above = x > thr;
ix = find(diff([0; above]) == 1);
tEdge = t_s(ix);

if isempty(tEdge), return; end
keep = true(size(tEdge));
for k = 2:numel(tEdge)
    if (tEdge(k) - tEdge(find(keep,1,'last'))) < minISI_s
        keep(k) = false;
    end
end
tEdge = tEdge(keep);
end

function [tEdgeO, thr] = detect_onebox_adc_edges(streamADC, adcIdx, nCh, Fs, opts)
datFile = fullfile(streamADC,'continuous.dat');
if exist(datFile,'file') ~= 2
    error('Missing %s', datFile);
end

% chunk size control
maxBytes = opts.maxChunkMB * 1024^2;
chunkN = floor(maxBytes / (2*double(nCh)));
chunkN = max(1, chunkN);
chunkN = min(chunkN, floor(opts.chunk_s * Fs));

% estimate threshold from first ~5 seconds
nInit = min(chunkN, floor(5*Fs));
fid = fopen(datFile,'r');
if fid < 0, error('Cannot open %s', datFile); end
X0 = fread(fid, [nCh nInit], 'int16=>double');
fclose(fid);
if isempty(X0), error('Empty OneBox ADC file'); end
x0 = X0(adcIdx,:);

thr = opts.thr_onebox;
if isempty(thr)
    thr = (prctile(x0,5) + prctile(x0,95))/2;
end

% scan file for rising edges
info = dir(datFile);
nSamples = floor(double(info.bytes) / (2*double(nCh)));

fid = fopen(datFile,'r');
if fid < 0, error('Cannot open %s', datFile); end

prevAbove = false;
lastEdgeSamp = -Inf;
minISI_samp = round(opts.minISI_s * Fs);

edges = zeros(1, 1e5); % prealloc, will grow if needed
nE = 0;

s0 = 1;
while s0 <= nSamples
    nRead = min(chunkN, nSamples - s0 + 1);
    X = fread(fid, [nCh nRead], 'int16=>double');
    if isempty(X), break; end
    x = X(adcIdx,:);

    above = (x > thr);
    d = diff([prevAbove, above]);
    ix = find(d == 1); % rising indices inside chunk (1..nRead)

    for k = 1:numel(ix)
        samp = (s0 - 1) + ix(k);
        if samp - lastEdgeSamp >= minISI_samp
            nE = nE + 1;
            if nE > numel(edges), edges = [edges, zeros(1, 1e5)]; end %#ok<AGROW>
            edges(nE) = samp;
            lastEdgeSamp = samp;
        end
    end

    prevAbove = above(end);
    s0 = s0 + nRead;
end

fclose(fid);

edges = edges(1:nE);
tEdgeO = (double(edges) - 1) / double(Fs); % seconds relative to segment start
end

function [tM2, tO2] = align_by_shift(tM, tO, sh)
if sh > 0
    tOa = tO(1+sh:end);
    tMa = tM;
elseif sh < 0
    tMa = tM(1-sh:end);
    tOa = tO;
else
    tMa = tM;
    tOa = tO;
end
n = min(numel(tMa), numel(tOa));
tM2 = tMa(1:n);
tO2 = tOa(1:n);
end

function streamRoot = oe_find_stream(segFolder, key)
streamRoot = '';
cand = dir(fullfile(segFolder,'recording*','continuous',['*' key '*']));
if isempty(cand)
    cand = dir(fullfile(segFolder,'Record Node*','experiment*','recording*','continuous',['*' key '*']));
end
for i = 1:numel(cand)
    p = fullfile(cand(i).folder, cand(i).name);
    if exist(fullfile(p,'continuous.dat'),'file')
        streamRoot = p;
        return
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

key = ['"stream_name": "' streamShort '"'];
ix = strfind(txt, key);
if isempty(ix)
    key = ['"stream_name":"', streamShort, '"'];
    ix = strfind(txt, key);
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