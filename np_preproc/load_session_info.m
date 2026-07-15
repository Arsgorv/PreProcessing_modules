function [fs, nCh, lsb, ycoords] = load_session_info(basepath)
%LOAD_SESSION_INFO  Sampling rate, channel count, bit->uV (and y-coords).
%   [fs,nCh,lsb,ycoords] = load_session_info(basepath)
%
%   Reads a CellExplorer 'continuous.session.mat' if present in basepath;
%   otherwise falls back to the OpenEphys 'structure.oebin' for the continuous
%   stream in basepath (searched up to 4 folders up). Lets the pipeline run on
%   sessions that were never processed by CellExplorer.
%
%   basepath : folder containing continuous.dat
%   ycoords  : per-channel depth (um) if available in session.mat, else []

ycoords = [];
sf = fullfile(basepath,'continuous.session.mat');
if exist(sf,'file')
    S = load(sf,'session'); e = S.session.extracellular;
    fs = e.srLfp; nCh = e.nChannels; lsb = e.leastSignificantBit;
    if isfield(e,'chanCoords') && isfield(e.chanCoords,'y'), ycoords = e.chanCoords.y(:); end
    return
end

% ---- fall back to OpenEphys structure.oebin ----
oebin = ''; p = basepath;
for up = 1:4
    cand = fullfile(p,'structure.oebin');
    if exist(cand,'file'), oebin = cand; break; end
    p = fileparts(p);
end
if isempty(oebin)
    error('load_session_info:notfound', ...
        'No continuous.session.mat and no structure.oebin found near\n  %s', basepath);
end

j = jsondecode(fileread(oebin));
cont = j.continuous; if ~iscell(cont), cont = num2cell(cont); end

% Pick the right continuous stream. OneBox/Neuropixels record an AP (~30 kHz)
% and an LFP (~2.5 kHz) stream per probe; the folder we are reading tells us
% which we want. Restrict to streams of THIS probe (folder_name minus the
% -AP/-LFP suffix), then choose AP = highest / LFP = lowest sample rate. This
% is robust even when the oebin folder_name differs from the on-disk folder.
[~,streamFolder] = fileparts(basepath); sf = lower(streamFolder);
probeTok = regexprep(sf,'[-_.](ap|lfp)$','');        % e.g. 'onebox-129.probea'
srs = zeros(1,numel(cont)); fns = cell(1,numel(cont));
for i = 1:numel(cont)
    srs(i) = double(cont{i}.sample_rate);
    fn = ''; if isfield(cont{i},'folder_name'), fn = cont{i}.folder_name; end
    fns{i} = lower(regexprep(fn,'[\\/]+$',''));
end
cand = find(~cellfun(@isempty, strfind(fns, probeTok)));
if isempty(cand) || isempty(probeTok), cand = 1:numel(cont); end
if ~isempty(strfind(sf,'lfp')),     [~,j] = min(srs(cand));
elseif ~isempty(strfind(sf,'ap')),  [~,j] = max(srs(cand));
else,                               j = 1;
end
sel = cont{cand(j)};

fs  = double(sel.sample_rate);
nCh = double(sel.num_channels);
bv  = arrayfun(@(c) double(c.bit_volts), sel.channels);
lsb = bv(1);

% sanity: num_channels must divide continuous.dat evenly, else wrong stream
df = fullfile(basepath,'continuous.dat');
if exist(df,'file')
    b = dir(df);
    if mod(b.bytes/2, nCh) ~= 0
        warning('load_session_info:nCh', ...
            'num_channels (%d) does not divide continuous.dat evenly - stream/oebin mismatch?', nCh);
    end
end
end
