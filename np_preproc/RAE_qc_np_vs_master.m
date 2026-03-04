function QC = RAE_qc_np_vs_master(datapath, opts)
% Compare per-segment duration in master (manifest) vs NP ProbeA-LFP.
% Optionally estimate OneBox-vs-master drift from 1Hz TTL recorded on AB channel (e.g. LFP91).

if nargin < 2, opts = struct(); end
if ~isfield(opts,'ttl_lfp_chan'), opts.ttl_lfp_chan = 91; end

TsRate = 1e4;

S = load(fullfile(datapath,'analysis','run_manifest_RAexp.mat'),'RUN');
seg = S.RUN.segments;

% load TTL LFP if present (stitched master axis)
ttlFile = fullfile(datapath,'ephys','LFPData', sprintf('LFP%d.mat', opts.ttl_lfp_chan));
haveTTL = exist(ttlFile,'file')==2;
if haveTTL
    Sttl = load(ttlFile,'LFP');
    Lttl = Sttl.LFP; % tsd
end

rows = {};
for i = 1:numel(seg)
    segFolder = fullfile(datapath,'ephys', seg(i).name);
    streamNP = oe_find_stream(segFolder, 'ProbeA-LFP');
    if isempty(streamNP)
        rows(end+1,:) = {seg(i).name, NaN, NaN, NaN, NaN, NaN}; %#ok<AGROW>
        continue
    end

    [FsNP, nChNP] = oe_read_stream_info(streamNP, NaN, NaN);
    if ~isfinite(FsNP), FsNP = 2500; end
    if ~isfinite(nChNP), nChNP = 384; end

    nSampNP = oe_nSamples_stream(streamNP, nChNP); % robust (timestamps.npy if exists, else bytes)
    durNP_s = double(nSampNP) / double(FsNP);

    durMaster_s = seg(i).t1_s - seg(i).t0_s;
    delta_s = durNP_s - durMaster_s;

    % predicted overlap at 1250 Hz
    overlap_samp_1250 = round(delta_s * 1250);

    % TTL drift estimate inside this segment (master axis)
    dtPulse = NaN;
    if haveTTL
        t0_ts = seg(i).t0_s * TsRate;
        t1_ts = seg(i).t1_s * TsRate;
        Lw = Restrict(Lttl, intervalSet(t0_ts, t1_ts));
        x = Data(Lw);
        t = Range(Lw) / TsRate;

        if ~isempty(x)
            thr = median(x) + 5*mad(x,1); % robust
            above = (x > thr);
            edge = find(diff([0; above])==1); % rising edges
            if numel(edge) > 5
                tt = t(edge);
                dt = diff(tt);
                dtPulse = median(dt); % should be ~1.0 s (OneBox seconds) measured in master seconds
            end
        end
    end

    rows(end+1,:) = {seg(i).name, durMaster_s, durNP_s, delta_s*1000, overlap_samp_1250, dtPulse}; 
end

QC = cell2table(rows, 'VariableNames', {'segment','dur_master_s','dur_np_s','delta_ms','overlap_samp_1250','ttl_median_dt_s'});
disp(QC)
end

% --- helpers (copy from your file or put in a shared utils file) ---
function n = oe_nSamples_stream(streamRoot, nCh)
% Prefer timestamps.npy length (fast header read), else use bytes.
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
% Robust to stream folders containing dots: "OneBox-102.ProbeA-LFP"
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

% stream folder name with dots
streamFolder = regexp(streamRoot, '[^\\/]+$', 'match', 'once');
if isempty(streamFolder), return; end

dotix = find(streamFolder=='.', 1, 'last');
if ~isempty(dotix)
    streamShort = streamFolder(dotix+1:end); % e.g. ProbeA-LFP
else
    streamShort = streamFolder;
end

key = ['"stream_name": "' streamShort '"'];
ix = strfind(txt, key);
if isempty(ix)
    key = ['"stream_name":"', streamShort, '"'];
    ix = strfind(txt, key);
    if isempty(ix), return; end
end
ix = ix(1);

w0 = max(1, ix-6000);
w1 = min(numel(txt), ix+6000);
win = txt(w0:w1);
pos = ix - w0 + 1;

% sample_rate nearest BEFORE stream_name
[s0, ~, tok] = regexp(win, 'sample_rate\"?\s*:\s*([0-9]+\.?[0-9]*)', 'start','end','tokens');
if ~isempty(s0)
    k = find(s0 < pos, 1, 'last');
    if isempty(k), k = 1; end
    Fs = str2double(tok{k}{1});
end

% num_channels nearest AFTER stream_name (else before)
[s1, ~, tok] = regexp(win, 'num_channels\"?\s*:\s*([0-9]+)', 'start','end','tokens');
if ~isempty(s1)
    k = find(s1 > pos, 1, 'first');
    if isempty(k), k = find(s1 < pos, 1, 'last'); end
    if isempty(k), k = 1; end
    nCh = str2double(tok{k}{1});
end
end
