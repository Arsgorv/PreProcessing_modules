function RUN = RAE_make_run_manifest(datapath)
% RAE_make_run_manifest
% Builds high-level run intervals for RA experiment sessions:
%   PreSleep, Conditioning, PostSleep, PostTest (+ optional _1/_2 segments)
% Assumes your stitched LFPData concatenates OE recordings in chronological order
% (no inserted gaps). Provides run boundaries in seconds and ticks (1e4).
%
% Saves:
%   <datapath>/analysis/run_manifest_RAexp.mat

TsRate = 1e4;

ephysDir = fullfile(datapath,'ephys');
if ~exist(ephysDir,'dir')
    error('RAE_make_run_manifest:NoEphys','No ephys folder: %s', ephysDir);
end

% find candidate phase folders
D = dir(ephysDir);
names = {D([D.isdir]).name};
names = names(~ismember(names,{'.','..','LFPData'}));

% keep only folders containing _RA_<Phase>
phList = {'PreSleep','Conditioning','PostSleep','PostTest'};
seg = struct('name',{},'phase',{},'segIdx',{},'tStartName',{},'path',{},'dur_s',{});

for i = 1:numel(names)
    nm = names{i};
    ph = '';
    for p = 1:numel(phList)
        if ~isempty(regexp(nm, ['_RA_' phList{p} '(\_\d+)?$'], 'once'))
            ph = phList{p};
            break
        end
    end
    if isempty(ph), continue; end

    % segment index suffix (_1/_2/...) default 0
    tok = regexp(nm, ['_RA_' ph '_(\d+)$'], 'tokens','once');
    if isempty(tok), segIdx = 0; else, segIdx = str2double(tok{1}); end

    % try parse timestamp in name: yyyy-mm-dd_HH-MM-SS
    tStartName = NaN;
    tok2 = regexp(nm, '(\d{4}-\d{2}-\d{2})_(\d{2}-\d{2}-\d{2})', 'tokens','once');
    if ~isempty(tok2)
        try
            tStartName = datenum([tok2{1} ' ' strrep(tok2{2},'-',':')], 'yyyy-mm-dd HH:MM:SS');
        catch
        end
    end

    seg(end+1).name = nm; 
    seg(end).phase = ph;
    seg(end).segIdx = segIdx;
    seg(end).tStartName = tStartName;
    seg(end).path = fullfile(ephysDir,nm);
    seg(end).dur_s = oe_segment_duration_seconds(seg(end).path);
end

if isempty(seg)
    error('RAE_make_run_manifest:NotRAExp','No *_RA_<Phase> folders found in %s', ephysDir);
end

% sort chronologically: prefer parsed timestamp, fallback to folder datenum
tSort = nan(numel(seg),1);
for i = 1:numel(seg)
    if isfinite(seg(i).tStartName)
        tSort(i) = seg(i).tStartName;
    else
        d = dir(seg(i).path);
        tSort(i) = d(1).datenum;
    end
end
[~,ord] = sort(tSort);
seg = seg(ord);

% accumulate offsets
t0 = 0;
for i = 1:numel(seg)
    seg(i).t0_s = t0;
    seg(i).t1_s = t0 + seg(i).dur_s;
    t0 = seg(i).t1_s;
end

% build phase-level runs (union of segments of same phase)
RUN = struct();
RUN.datapath = datapath;
RUN.TsRate = TsRate;
RUN.segments = seg;

for p = 1:numel(phList)
    ph = phList{p};
    ii = find(strcmp({seg.phase}, ph));
    if isempty(ii)
        RUN.run.(ph).t0_s = NaN;
        RUN.run.(ph).t1_s = NaN;
        continue
    end
    RUN.run.(ph).t0_s = seg(ii(1)).t0_s;
    RUN.run.(ph).t1_s = seg(ii(end)).t1_s;
    RUN.run.(ph).t0_ts = RUN.run.(ph).t0_s * TsRate;
    RUN.run.(ph).t1_ts = RUN.run.(ph).t1_s * TsRate;
end

% save
outDir = fullfile(datapath,'analysis');
if ~exist(outDir,'dir'), mkdir(outDir); end
save(fullfile(outDir,'run_manifest_RAexp.mat'), 'RUN');

end

% ---------------- helpers (no nesting, local subfunctions are OK) ----------------
function dur_s = oe_segment_duration_seconds(segPath)
% robust to both pre/post fix_folder_structure layouts

% try "fixed" layout first:
cand = dir(fullfile(segPath,'recording*','continuous','Acquisition_Board*','timestamps.npy'));
if isempty(cand)
    % try "unfixed" layout:
    cand = dir(fullfile(segPath,'Record Node*','experiment*','recording*','continuous','Acquisition_Board*','timestamps.npy'));
end

Fs = 30000; % fallback
Fs2 = try_read_sample_rate(segPath);
if isfinite(Fs2) && Fs2 > 1, Fs = Fs2; end

if ~isempty(cand)
    tsPath = fullfile(cand(1).folder, cand(1).name);
    n = npy_numel(tsPath);
    dur_s = double(n) / double(Fs);
    return
end

% fallback: use continuous.dat size with channel count if timestamps.npy missing
cand2 = dir(fullfile(segPath,'recording*','continuous','Acquisition_Board*','continuous.dat'));
if isempty(cand2)
    cand2 = dir(fullfile(segPath,'Record Node*','experiment*','recording*','continuous','Acquisition_Board*','continuous.dat'));
end
if isempty(cand2)
    error('oe_segment_duration_seconds:NoOE','Cannot find timestamps.npy or continuous.dat in %s', segPath);
end

nChan = try_read_num_channels(segPath);
if ~isfinite(nChan) || nChan < 1
    nChan = 80; % your typical OB+FMA
end

bytes = cand2(1).bytes;
nSamples = floor(double(bytes) / (2*double(nChan))); % int16
dur_s = nSamples / double(Fs);

end

function n = npy_numel(npyPath)
fid = fopen(npyPath,'r');
if fid < 0, error('npy_numel:OpenFail','Cannot open %s', npyPath); end

magic = fread(fid,6,'uint8=>char')';
if ~strcmp(magic, char([147 'NUMPY']))
    fclose(fid); error('npy_numel:BadMagic','Not a npy file: %s', npyPath);
end
ver = fread(fid,2,'uint8');
if ver(1) == 1
    hlen = fread(fid,1,'uint16');
else
    hlen = fread(fid,1,'uint32');
end
hdr = fread(fid, double(hlen), 'uint8=>char')';
fclose(fid);

% parse shape from header string, e.g. "(12345,)"
tok = regexp(hdr, 'shape\W*\(\s*([0-9]+)', 'tokens','once');
if isempty(tok)
    error('npy_numel:NoShape','Cannot parse shape in header for %s', npyPath);
end
n = str2double(tok{1});
end

function Fs = try_read_sample_rate(segPath)
Fs = NaN;
% look for structure.oebin
cand = dir(fullfile(segPath,'recording*','structure.oebin'));
if isempty(cand)
    cand = dir(fullfile(segPath,'Record Node*','experiment*','recording*','structure.oebin'));
end
if isempty(cand), return; end

txt = '';
try
    txt = fileread(fullfile(cand(1).folder, cand(1).name));
catch
    return
end

% regexp-based (robust to format differences)
% take first occurrence near "Acquisition_Board"
ix = strfind(txt, 'Acquisition_Board');
if isempty(ix)
    ix = strfind(txt, 'Acquisition Board');
end
if isempty(ix)
    ix = 1;
else
    ix = ix(1);
end

win = txt(ix:min(ix+5000, numel(txt)));
tok = regexp(win, 'sample_rate\"?\s*:\s*([0-9]+\.?[0-9]*)', 'tokens','once');
if isempty(tok), return; end
Fs = str2double(tok{1});
end

function nChan = try_read_num_channels(segPath)
nChan = NaN;
cand = dir(fullfile(segPath,'recording*','structure.oebin'));
if isempty(cand)
    cand = dir(fullfile(segPath,'Record Node*','experiment*','recording*','structure.oebin'));
end
if isempty(cand), return; end

txt = '';
try
    txt = fileread(fullfile(cand(1).folder, cand(1).name));
catch
    return
end

ix = strfind(txt, 'Acquisition_Board');
if isempty(ix)
    ix = strfind(txt, 'Acquisition Board');
end
if isempty(ix)
    ix = 1;
else
    ix = ix(1);
end

win = txt(ix:min(ix+5000, numel(txt)));
tok = regexp(win, 'num_channels\"?\s*:\s*([0-9]+)', 'tokens','once');
if isempty(tok), return; end
nChan = str2double(tok{1});
end