function Master_LFP_NP_preproc(sessions, opts)
% Master_LFP_NP_preproc
% Export Neuropixels LFP (ProbeA-LFP) onto the MASTER (Acquisition Board) stitched timeline.
%
% Outputs:
%   per_channel:
%     <datapath>/ephys/LFPDataNP/LFP<ch>.mat      (variable LFP as tsd)
%   matrix (recommended for many channels):
%     <datapath>/ephys/LFPDataNP/LFPmat_<seg>.mat (variables Y int16, t0_ts, dt_ts, nOut, chUse, Fs_ds)
%
% Required:
%   opts.np_channels : 1-based channels within ProbeA-LFP stream
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

TsRate = 1e4;

for s = 1:numel(sessions)
    datapath = sessions{s};
    disp('------------------------------------------')
    disp(['[Master_LFP_NP_preproc] ' datapath])

    outDir = fullfile(datapath,'ephys','LFPDataNP');
    if ~exist(outDir,'dir'), mkdir(outDir); end

    % ---------------- segments on MASTER axis ----------------
    mf = fullfile(datapath,'analysis','run_manifest_RAexp.mat');
    if exist(mf,'file') == 2
        Srun = load(mf,'RUN'); RUN = Srun.RUN;
        segs = RUN.segments;
        for i = 1:numel(segs)
            segs(i).t0_ts = segs(i).t0_s * TsRate;
            segs(i).t1_ts = segs(i).t1_s * TsRate;
        end
    else
        error('Master_LFP_NP_preproc:NeedManifest','Missing run_manifest_RAexp.mat. Required for correct master boundaries.');
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
    matFiles = {};

    % =====================================================================
    % MATRIX MODE (stream-write to avoid huge RAM; safe for 384 channels)
    % =====================================================================
    if strcmp(save_style_use,'matrix')

        for i = 1:numel(segs)
            segFolder = fullfile(datapath,'ephys', segs(i).name);
            streamRoot = oe_find_stream(segFolder, 'ProbeA-LFP');
            if isempty(streamRoot)
                warning('  no ProbeA-LFP in %s', segs(i).name);
                continue
            end

            [Fs, nCh] = oe_read_stream_info(streamRoot, NaN, NaN);
            if ~isfinite(Fs), Fs = 2500; end
            if ~isfinite(nCh), nCh = 384; end
            Rraw = Fs / opts.lfp_fs;
            R = round(Rraw);
            if abs(Rraw - R) > 1e-6
                warning('Fs=%.3f not integer divisible by targetFs=%.3f. Using R=%d (Fs_ds=%.3f).', Fs, opts.lfp_fs, R, Fs/R);
            end

            keep = (chList >= 1) & (chList <= nCh);
            if ~any(keep)
                warning('  none of requested channels exist in %s', segs(i).name);
                continue
            end
            chUse = chList(keep);

            outMat = fullfile(outDir, ['LFPmat_' segs(i).name '.mat']);
            if exist(outMat,'file')==2 && ~opts.force
                disp(['  skip existing mat: ' outMat]);
                matFiles{end+1} = outMat; %#ok<AGROW>
                continue
            end

            segStart_ts = double(segs(i).t0_ts);
            segStop_ts  = double(segs(i).t1_ts);

            [Fs_ds, dt_ts, nOut] = oe_write_lfp_matrix_fread( ...
                streamRoot, outMat, chUse, nCh, Fs, opts.lfp_fs, ...
                opts.maxChunkMB, opts.chunk_s, segStart_ts, segStop_ts, TsRate);

            t0_ts = segStart_ts; 
            save(outMat, 't0_ts','dt_ts','nOut','chUse','Fs_ds','-append');

            Fs_ds_bySeg(i) = Fs_ds;
            matFiles{end+1} = outMat; %#ok<AGROW>
            disp(['  saved mat: ' outMat]);
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
            streamRoot = oe_find_stream(segFolder, 'ProbeA-LFP');
            if isempty(streamRoot)
                warning('  no ProbeA-LFP in %s', segs(i).name);
                continue
            end

            [Fs, nCh] = oe_read_stream_info(streamRoot, NaN, NaN);
            if ~isfinite(Fs), Fs = 2500; end
            if ~isfinite(nCh), nCh = 384; end
            
            Rraw = Fs / opts.lfp_fs;
            R = round(Rraw);
            if abs(Rraw - R) > 1e-6
                warning('Fs=%.3f not integer divisible by targetFs=%.3f. Using R=%d (Fs_ds=%.3f).', Fs, opts.lfp_fs, R, Fs/R);
            end

            keep = (chList >= 1) & (chList <= nCh);
            if ~any(keep), continue; end
            chUse = chList(keep);

            segStart_ts = double(segs(i).t0_ts);
            segStop_ts  = double(segs(i).t1_ts);

            % integer downsample
            R = max(1, round(Fs / opts.lfp_fs));
            Fs_ds = Fs / R;

            dt_ts_real = TsRate / Fs_ds;
            dt_ts = round(dt_ts_real);
            if abs(dt_ts_real - dt_ts) > 1e-9
                warning('dt_ts not integer (TsRate/Fs_ds). TsRate=%g Fs_ds=%g dt_ts=%g', TsRate, Fs_ds, dt_ts_real);
            end

            % max samples allowed so that t(end) < segStop_ts
            maxNOut = floor((segStop_ts - segStart_ts) / dt_ts);
            if maxNOut < 1
                warning('Bad segment boundaries: %s (stop<=start)', segs(i).name);
                continue
            end

            [Yseg, Fs_ds2] = oe_extract_channels_lfp_fread( ...
                streamRoot, chUse, nCh, Fs, opts.lfp_fs, opts.maxChunkMB, opts.chunk_s, maxNOut);

            if isempty(Yseg), continue; end
            Fs_ds_bySeg(i) = Fs_ds2;

            % build time from MASTER start
            t_ts = segStart_ts + (0:size(Yseg,1)-1)' * dt_ts;

            % append parts per channel, trimming any accidental overlap
            for jj = 1:numel(chUse)
                ch = chUse(jj);
                j = find(chList == ch, 1, 'first');

                y = Yseg(:,jj);

                if ~isempty(parts{j})
                    prev = Range(parts{j}{end});
                    prevEnd = prev(end);

                    if t_ts(1) <= prevEnd
                        kk = (t_ts > prevEnd);
                        t_ts2 = t_ts(kk);
                        y2 = y(kk);
                    else
                        t_ts2 = t_ts;
                        y2 = y;
                    end
                else
                    t_ts2 = t_ts;
                    y2 = y;
                end

                if ~isempty(t_ts2)
                    parts{j}{end+1} = tsd(t_ts2, y2);
                end
            end
        end

        % save per channel
        for j = 1:nSel
            ch = chList(j);
            outFile = fullfile(outDir, ['LFP' num2str(ch) '.mat']);
            if exist(outFile,'file')==2 && ~opts.force
                disp(['  skip existing: ' outFile]);
                continue
            end
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
            save(outFile,'LFP','-v7.3');
            disp(['  saved: ' outFile]);
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

    save(fullfile(outDir,'InfoLFP_NP.mat'), 'Info');

end

end

% ============================ helpers ============================

function segs = list_ephys_segments_simple(datapath)
ephysDir = fullfile(datapath,'ephys');
D = dir(ephysDir);
D = D([D.isdir]);
D = D(~ismember({D.name},{'.','..','LFPData','LFPDataNP'}));
[~,ord] = sort([D.datenum]);
D = D(ord);
segs = struct('name',{},'t0_s',{},'t0_ts',{},'dur_s',{},'t1_s',{},'t1_ts',{});
for i = 1:numel(D)
    segs(i).name = D(i).name; 
end
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

streamFolder = regexp(streamRoot, '[^\\/]+$', 'match', 'once');
if isempty(streamFolder), return; end

dotix = find(streamFolder=='.', 1, 'last');
if ~isempty(dotix)
    streamShort = streamFolder(dotix+1:end);
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
% Fast chunked extraction using fread on contiguous blocks.
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
if nSel==0
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

maxBytes = maxChunkMB * 1024^2;
chunkN = floor(maxBytes / (2*double(nCh)));
chunkN = max(R, floor(chunkN/R)*R);
chunkN = max(chunkN, Fs * 2);
chunkN = min(chunkN, Fs * 120);

if isfinite(chunk_s) && chunk_s > 0
    chunkN2 = floor(chunk_s * Fs);
    chunkN2 = max(R, floor(chunkN2/R)*R);
    chunkN = min(chunkN, chunkN2);
end

nOutEst = floor(nSamples / R);
Y_ds = zeros(nOutEst, nSel, 'int16');

fid = fopen(datFile,'r');
if fid < 0
    Y_ds = []; Fs_ds = NaN; return
end

pos = 1;
s0 = 1;
while s0 <= nSamples
    nRead = min(chunkN, nSamples - s0 + 1);

    X = fread(fid, [nCh nRead], 'int16=>int16');
    if isempty(X), break; end

    Xs = X(chList, :);
    Xd = Xs(:, 1:R:end);
    k  = size(Xd,2);

    if k > 0
        Y_ds(pos:pos+k-1, :) = Xd.';
        pos = pos + k;
    end

    s0 = s0 + nRead;
end

fclose(fid);
Y_ds = Y_ds(1:pos-1, :);

if isfinite(maxNOut) && maxNOut > 0 && size(Y_ds,1) > maxNOut
    Y_ds = Y_ds(1:maxNOut, :);
end
end

function [Fs_ds, dt_ts, nOut] = oe_write_lfp_matrix_fread(streamRoot, outMat, chList, nCh, Fs, targetFs, maxChunkMB, chunk_s, segStart_ts, segStop_ts, TsRate)
% Stream-write Y into outMat, clamped to master segment stop.

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
nOutEst = floor(nSamplesLimit / R);

nSel = numel(chList);

M = matfile(outMat, 'Writable', true);
M.Y = zeros(nOutEst, nSel, 'int16');

maxBytes = maxChunkMB * 1024^2;
chunkN = floor(maxBytes / (2*double(nCh)));
chunkN = max(R, floor(chunkN/R)*R);
chunkN = max(chunkN, Fs * 2);
chunkN = min(chunkN, Fs * 120);

if nargin < 8 || isempty(chunk_s), chunk_s = 60; end
if isfinite(chunk_s) && chunk_s > 0
    chunkN2 = floor(chunk_s * Fs);
    chunkN2 = max(R, floor(chunkN2/R)*R);
    chunkN = min(chunkN, chunkN2);
end

fid = fopen(datFile,'r');
if fid < 0, error('oe_write_lfp_matrix_fread:OpenFail','Cannot open %s', datFile); end

pos = 1;
s0 = 1;
while s0 <= nSamplesLimit
    nRead = min(chunkN, nSamplesLimit - s0 + 1);

    X = fread(fid, [nCh nRead], 'int16=>int16');
    if isempty(X), break; end

    Xs = X(chList, :);
    Xd = Xs(:, 1:R:end);
    k  = size(Xd,2);

    if k > 0
        M.Y(pos:pos+k-1, :) = Xd.';
        pos = pos + k;
    end

    s0 = s0 + nRead;
end

fclose(fid);

nOut = pos - 1;
end