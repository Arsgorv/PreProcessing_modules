ExtractSessionMetadata_OB({'Shropshire'})
PlotSessionMetadata_OB(sessions, save_dir)

function ExtractSessionMetadata_OB(animals, setups)
% ExtractSessionMetadata_OB  Build per-session metadata for OB ferret recordings.
%
% GOAL
%   For each session of Brynza / Labneh / Shropshire (freely-moving and
%   head-fixed), collect a metadata structure and save it as
%   <session_folder>/session_metadata.mat. The structure contains:
%       animal              - ferret name (derived from the data path)
%       animal_age_days     - age at recording (recording date - birth date)
%       animal_weight_g     - body weight interpolated to the recording date
%       mean_temperature_c  - mean stabulation temp during the recording
%       recording_date      - 'yyyy-mm-dd'
%       start_time          - datetime (from OpenEphys settings.xml)
%       end_time            - datetime (start_time + duration_sec)
%       duration_sec        - recorded data length (from sample count)
%       recording_path      - full path to the session folder
%   plus QC / provenance fields (time_source, sampling_rate_hz, etc.).
%
% PROVENANCE
%   - sessions          : PathForExperimentsOB (must be on the MATLAB path).
%   - recording_date    : <DATE> tag in settings.xml (e.g. "23 Jan 2024
%                         18:27:33"). Fallbacks: ExpeInfo.date, folder name,
%                         .dat filename.
%   - start_time        : same <DATE> tag (wall-clock, second-resolution).
%                         Fallbacks: TimeRec.mat (TimeBeginRec/Allfiles),
%                         behavResources.mat, structure.oebin mtime.
%   - duration_sec      : sum over recording*/structure.oebin of
%                         (continuous.dat bytes)/(2*num_channels)/sample_rate.
%                         Gold standard for recorded data length.
%   - end_time          : start_time + duration_sec.
%   - animal_weight_g   : Z:\Arsenii\ferretware\poids_biologid.csv
%                         (';'-delimited, no header). Interpolated to the
%                         recording date within +/- 7 days (linear; nearest
%                         within window if not bracketed; mean if same-day).
%   - birth date / age  : Z:\Arsenii\ferretware\Suivi information furets.xlsx
%                         Auto-detected name column ("Nom") and birth-date
%                         column ("Date de naissance", French strings like
%                         "3 juin 2020"). Falls back to col 1 / col 7.
%   - mean_temperature_c: \\129.199.81.18\data5\Arsenii\ferretware\Temperature.txt
%                         tab-delimited; mean over [start_time, end_time]
%                         (padded by +/- 1 h if window is empty).
%
% USAGE
%   ExtractSessionMetadata_OB
%   ExtractSessionMetadata_OB({'Brynza'})
%   ExtractSessionMetadata_OB({'Labneh'}, {'freely-moving'})
%
% MATLAB R2018b conventions (no string/contains/readcell/readtable detect).
% Arsenii / OB project.

%% -------------------- parameters / defaults --------------------
if nargin < 1 || isempty(animals)
    animals = {'Brynza', 'Labneh', 'Shropshire'};
end
if nargin < 2 || isempty(setups)
    setups = {'freely-moving', 'head-fixed'};
end
if ~iscell(animals), animals = {animals}; end
if ~iscell(setups),  setups  = {setups};  end

weight_csv  = 'Z:\Arsenii\ferretware\poids_biologid.csv';
birth_xlsx  = 'Z:\Arsenii\ferretware\Suivi information furets.xlsx';
temp_txt    = '\\129.199.81.18\data5\Arsenii\ferretware\Temperature.txt';

weight_window_days = 7;       % +/- window for weight interpolation
temp_pad_h         = 1;       % +/- hour padding when temp samples are sparse
dur_mismatch_warn  = 120;     % warn if |duration_sec - TimeRec span| > this (s)

%% -------------------- build the session list --------------------
sess_paths = {};
sess_expe  = {};
for s = 1:numel(setups)
    D = PathForExperimentsOB(animals, setups{s});
    sess_paths = [sess_paths, D.path];    %#ok<AGROW>
    sess_expe  = [sess_expe,  D.ExpeInfo];%#ok<AGROW>
end
[sess_paths, keep] = unique(sess_paths, 'stable');
sess_expe = sess_expe(keep);
nSess = numel(sess_paths);
fprintf('Found %d sessions for {%s} x {%s}.\n', nSess, ...
    strjoin(animals, ', '), strjoin(setups, ', '));

%% -------------------- load reference tables once --------------------
[wname, wdatenum, wval] = load_weight_table(weight_csv);
birthMap = load_birth_table(birth_xlsx);
[tdn, tval] = load_temperature_table(temp_txt);

%% -------------------- per-session loop --------------------
summary = cell(nSess, 9);
for i = 1:nSess
    p = sess_paths{i};
    E = sess_expe{i};
    fprintf('\n[%d/%d] %s\n', i, nSess, p);

    if exist(p, 'dir') ~= 7
        warning('Session folder not found - skipped: %s', p);
        continue;
    end

    animal = derive_animal(p, animals);
    if isempty(animal)
        warning('Could not derive animal name from path - skipped: %s', p);
        continue;
    end

    % --- date + start time of day (primary: settings.xml) ---
    [recDateNum, startHMS, timeSrc] = get_date_and_start(p, E);
    if isnan(recDateNum)
        warning('Could not determine recording date - skipped: %s', p);
        continue;
    end
    dateStr = datestr(recDateNum, 'yyyy-mm-dd');

    if all(~isnan(startHMS))
        startNum = recDateNum + (startHMS(1)*3600 + startHMS(2)*60 + startHMS(3)) / 86400;
    else
        startNum = NaN;
    end

    % --- duration from actual sample count (gold standard) ---
    [dur_sec, srUsed] = get_duration_from_data(p);

    % --- end time = start + duration ---
    if ~isnan(startNum) && ~isnan(dur_sec)
        endNum = startNum + dur_sec / 86400;
    else
        endNum = NaN;
    end

    % --- QC: compare with TimeRec.mat span if it exists ---
    [tr_start, tr_end] = read_timerec(p);
    if ~isnan(tr_start) && ~isnan(tr_end) && ~isnan(dur_sec)
        tr_span = (tr_end - tr_start) * 86400;
        if abs(tr_span - dur_sec) > dur_mismatch_warn
            fprintf(['  [QC] sample duration (%.1f s) and TimeRec span (%.1f s) ' ...
                     'differ by %.1f s.\n'], dur_sec, tr_span, abs(tr_span - dur_sec));
        end
    end

    % --- weight, interpolated ---
    [weight_g, nUsed] = interp_weight(wname, wdatenum, wval, animal, recDateNum, weight_window_days);

    % --- mean temperature during recording window ---
    temp_c = mean_temperature(tdn, tval, startNum, endNum, temp_pad_h);

    % --- age in days ---
    if isKey(birthMap, lower(animal))
        birthNum = birthMap(lower(animal));
        age_days = recDateNum - birthNum;
    else
        birthNum = NaN; age_days = NaN;
        warning('No birth date found for %s in %s', animal, birth_xlsx);
    end

    % --- assemble structure ---
    session_metadata = struct();
    session_metadata.animal             = animal;
    session_metadata.animal_age_days    = age_days;
    session_metadata.animal_weight_g    = weight_g;
    session_metadata.mean_temperature_c = temp_c;
    session_metadata.recording_date     = dateStr;
    if ~isnan(startNum)
        session_metadata.start_time = datetime(startNum, 'ConvertFrom', 'datenum');
    else
        session_metadata.start_time = NaT;
    end
    if ~isnan(endNum)
        session_metadata.end_time = datetime(endNum, 'ConvertFrom', 'datenum');
    else
        session_metadata.end_time = NaT;
    end
    session_metadata.duration_sec   = dur_sec;
    session_metadata.recording_path = p;

    % --- QC / provenance ---
    session_metadata.duration_hms       = sec2hms(dur_sec);
    session_metadata.time_source        = timeSrc;
    session_metadata.sampling_rate_hz   = srUsed;
    if ~isnan(birthNum)
        session_metadata.birth_date = datestr(birthNum, 'yyyy-mm-dd');
    else
        session_metadata.birth_date = '';
    end
    session_metadata.weight_npoints_used = nUsed;

    outfile = fullfile(p, 'session_metadata.mat');
    save(outfile, 'session_metadata');
    fprintf('  saved: %s\n', outfile);
    fprintf('  animal=%s | date=%s | start=%s | dur=%s | weight=%.1fg | age=%dd | T=%.1f C\n', ...
        animal, dateStr, datestr_safe(session_metadata.start_time), ...
        session_metadata.duration_hms, weight_g, round_or_nan(age_days), temp_c);

    summary(i, :) = {animal, dateStr, datestr_safe(session_metadata.start_time), ...
        datestr_safe(session_metadata.end_time), session_metadata.duration_hms, ...
        weight_g, round_or_nan(age_days), timeSrc, p};
end

%% -------------------- console summary --------------------
fprintf('\n================ SUMMARY (%d sessions) ================\n', nSess);
fprintf('%-11s %-11s %-9s %-9s %-10s %8s %6s  %-14s\n', ...
    'animal', 'date', 'start', 'end', 'duration', 'weight', 'age', 'time_src');
for i = 1:nSess
    if isempty(summary{i, 1}), continue; end
    fprintf('%-11s %-11s %-9s %-9s %-10s %8.1f %6d  %-14s\n', ...
        summary{i,1}, summary{i,2}, summary{i,3}, summary{i,4}, ...
        summary{i,5}, summary{i,6}, summary{i,7}, summary{i,8});
end
fprintf('=======================================================\n');

end % main


%% ===================================================================
%% Local helper functions
%% ===================================================================

function animal = derive_animal(p, animals)
animal = '';
pl = lower(p);
for k = 1:numel(animals)
    if ~isempty(strfind(pl, lower(animals{k}))) %#ok<STREMP>
        animal = animals{k};
        return;
    end
end
end


function [recDateNum, startHMS, src] = get_date_and_start(p, E)
% PRIMARY: parse <DATE> from settings.xml at the session root (or recursively).
% FALLBACK: ExpeInfo.date + TimeRec.mat / behavResources.mat / OE file mtime.
recDateNum = NaN; startHMS = [NaN NaN NaN]; src = 'none';

% --- (1) settings.xml ----------------------------------------------------
cand = dir(fullfile(p, 'settings.xml'));
if isempty(cand)
    cand = dir(fullfile(p, '**', 'settings.xml'));
end
best = NaN;
for k = 1:numel(cand)
    fp = fullfile(cand(k).folder, cand(k).name);
    try
        txt = fileread(fp);
    catch
        continue;
    end
    tok = regexp(txt, '<DATE>\s*([^<]+?)\s*</DATE>', 'tokens', 'once');
    if isempty(tok), continue; end
    d = parse_settings_date(tok{1});
    if ~isnan(d) && (isnan(best) || d < best)
        best = d;
    end
end
if ~isnan(best)
    v = datevec(best);
    recDateNum = datenum(v(1), v(2), v(3));
    startHMS = v(4:6);
    src = 'settings.xml';
    return;
end

% --- (2) ExpeInfo.date ---------------------------------------------------
if isempty(E)
    eFile = fullfile(p, 'ExpeInfo.mat');
    if exist(eFile, 'file') == 2
        S = load(eFile, 'ExpeInfo');
        if isfield(S, 'ExpeInfo'), E = S.ExpeInfo; end
    end
end
if isstruct(E) && isfield(E, 'date') && ~isempty(E.date)
    d = E.date;
    if isnumeric(d), d = num2str(d); end
    d = strtrim(d);
    if length(d) == 6, d = ['20' d]; end
    if length(d) == 8
        recDateNum = datenum(d, 'yyyymmdd');
    end
end

% --- date fallback to folder name / .dat ---------------------------------
if isnan(recDateNum)
    [~, fname] = fileparts(p);
    tok = regexp(fname, '(\d{8})', 'tokens', 'once');
    if ~isempty(tok), recDateNum = datenum(tok{1}, 'yyyymmdd'); end
end
if isnan(recDateNum)
    dd = dir(fullfile(p, '*.dat'));
    for k = 1:numel(dd)
        tok = regexp(dd(k).name, '_(\d{8})_', 'tokens', 'once');
        if ~isempty(tok), recDateNum = datenum(tok{1}, 'yyyymmdd'); break; end
    end
end

% --- start HMS fallback: TimeRec.mat -> behavResources.mat -> OE mtime ---
[tr_start, ~] = read_timerec(p);
if ~isnan(tr_start)
    v = datevec(tr_start);
    if isnan(recDateNum)
        recDateNum = datenum(v(1), v(2), v(3));
    end
    startHMS = v(4:6);
    src = 'TimeRec.mat';
    return;
end

% behavResources.mat
f = fullfile(p, 'behavResources.mat');
if exist(f, 'file') == 2
    w = whos('-file', f);
    if any(strcmp({w.name}, 'TimeBeginRec'))
        S = load(f, 'TimeBeginRec');
        hms = extract_hms(S.TimeBeginRec(1, :));
        if all(~isnan(hms))
            startHMS = hms;
            src = 'behavResources.mat';
            return;
        end
    end
end

% structure.oebin mtime
oeb = dir(fullfile(p, '**', 'structure.oebin'));
if ~isempty(oeb)
    [~, ix] = min([oeb.datenum]);
    v = datevec(oeb(ix).datenum);
    if isnan(recDateNum)
        recDateNum = datenum(v(1), v(2), v(3));
    end
    startHMS = v(4:6);
    src = 'OE mtime';
end
end


function dn = parse_settings_date(s)
% Parse "23 Jan 2024 18:27:33" (OpenEphys settings.xml). Returns datenum or NaN.
dn = NaN;
s = strtrim(s);
tok = regexp(s, '(\d{1,2})\s+([A-Za-z]+)\s+(\d{4})\s+(\d{1,2}):(\d{1,2}):(\d{1,2})', 'tokens', 'once');
if isempty(tok), return; end
day_ = str2double(tok{1});
mon_ = parse_en_month(tok{2});
yr_  = str2double(tok{3});
H_   = str2double(tok{4});
M_   = str2double(tok{5});
S_   = str2double(tok{6});
if isnan(mon_) || any(isnan([day_ yr_ H_ M_ S_])), return; end
dn = datenum(yr_, mon_, day_, H_, M_, S_);
end


function m = parse_en_month(s)
months = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'};
ix = find(strncmpi(s, months, 3), 1);
if isempty(ix), m = NaN; else, m = ix; end
end


function [t_start, t_end] = read_timerec(p)
% Returns datenum of begin/end from TimeRec.mat (date from settings.xml/ExpeInfo
% applied externally). Here we only have [H M S], so we reconstruct on the
% recording's date externally - this returns absolute datenum using the date
% from the file modification time as a stand-in.
t_start = NaN; t_end = NaN;
f = fullfile(p, 'TimeRec.mat');
if exist(f, 'file') ~= 2, return; end
S = load(f);
if ~isfield(S, 'TimeBeginRec') || ~isfield(S, 'TimeEndRec'), return; end
hs = extract_hms(S.TimeBeginRec(1, :));
he = extract_hms(S.TimeEndRec(1, :));
if any(isnan(hs)) || any(isnan(he)), return; end

% Pull date from TimeBeginRec_Allfiles if it carries Y M D, else from the file
if isfield(S, 'TimeBeginRec_Allfiles') && size(S.TimeBeginRec_Allfiles, 2) >= 6
    vv = S.TimeBeginRec_Allfiles(1, :);
    dnum0 = datenum(vv(1), vv(2), vv(3));
else
    info = dir(f);
    v = datevec(info.datenum);
    dnum0 = datenum(v(1), v(2), v(3));
end
t_start = dnum0 + (hs(1)*3600 + hs(2)*60 + hs(3)) / 86400;
t_end   = dnum0 + (he(1)*3600 + he(2)*60 + he(3)) / 86400;
if t_end < t_start, t_end = t_end + 1; end
end


function [dur_sec, srUsed] = get_duration_from_data(p)
% Total recorded duration (s) from raw sample counts, summed across recordings.
dur_sec = 0; srUsed = NaN; found = false;
oeb = dir(fullfile(p, '**', 'structure.oebin'));
for k = 1:numel(oeb)
    oebPath = fullfile(oeb(k).folder, oeb(k).name);
    try
        info = jsondecode(fileread(oebPath));
    catch
        continue;
    end
    if ~isfield(info, 'continuous') || isempty(info.continuous), continue; end
    c = info.continuous(1);
    sr  = double(c.sample_rate);
    nch = double(c.num_channels);
    folder_name = strtrim(c.folder_name);
    datPath = fullfile(oeb(k).folder, 'continuous', folder_name, 'continuous.dat');
    dd = dir(datPath);
    if isempty(dd), continue; end
    nSamples = dd.bytes / (2 * nch);     % int16 = 2 bytes/sample/channel
    dur_sec  = dur_sec + nSamples / sr;
    srUsed   = sr;
    found    = true;
end
if ~found, dur_sec = NaN; end
end


function [w, nUsed] = interp_weight(names, wdatenum, wval, animal, recDateNum, win)
w = NaN; nUsed = 0;
mask = strcmpi(names, animal);
d = wdatenum(mask); v = wval(mask);
ok = ~isnan(d) & ~isnan(v); d = d(ok); v = v(ok);
if isempty(d), return; end

dday = floor(d);
[ud, ~, grp] = unique(dday);
uv = accumarray(grp, v, [], @mean);

sel = abs(ud - recDateNum) <= win;
ds = ud(sel); vs = uv(sel);
nUsed = numel(ds);
if isempty(ds), return; end
[ds, o] = sort(ds); vs = vs(o);

if any(ds == recDateNum)
    w = mean(vs(ds == recDateNum));
elseif recDateNum >= min(ds) && recDateNum <= max(ds)
    w = interp1(ds, vs, recDateNum, 'linear');
else
    [~, ix] = min(abs(ds - recDateNum));
    w = vs(ix);
end
end


function t = mean_temperature(tdn, tval, startNum, endNum, pad_h)
% Mean temperature inside [startNum, endNum]; extends by +/- pad_h if empty.
t = NaN;
if isempty(tdn) || isnan(startNum) || isnan(endNum), return; end
mask = tdn >= startNum & tdn <= endNum;
if ~any(mask) && pad_h > 0
    pad = pad_h / 24;
    mask = tdn >= (startNum - pad) & tdn <= (endNum + pad);
end
if any(mask), t = mean(tval(mask)); end
end


function [names, dnum, vals] = load_weight_table(csvfile)
% poids_biologid.csv: biologID;name;status;datetime;weight;unit  (no header).
names = {}; dnum = []; vals = [];
if exist(csvfile, 'file') ~= 2
    warning('Weight file not found: %s', csvfile);
    return;
end
fid = fopen(csvfile, 'r', 'n', 'UTF-8');
C = textscan(fid, '%s%s%s%s%f%s', 'Delimiter', ';');
fclose(fid);
names = C{2};
vals  = C{5};
dstr  = C{4};
% vectorized datenum
try
    dnum = datenum(dstr, 'yyyy-mm-dd HH:MM:SS');
catch
    dnum = NaN(numel(dstr), 1);
    for i = 1:numel(dstr)
        try, dnum(i) = datenum(dstr{i}, 'yyyy-mm-dd HH:MM:SS'); catch, dnum(i) = NaN; end
    end
end
fprintf('Loaded %d weight measurements from %s\n', numel(names), csvfile);
end


function birthMap = load_birth_table(xlsxfile)
% 'Suivi information furets.xlsx': autodetect "Nom" and "Date de naissance"
% column indices from the header row (so hidden columns B/D don't matter).
% Falls back to (col 1, col 7) -- matches the current sheet (column G).
birthMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
if exist(xlsxfile, 'file') ~= 2
    warning('Birth-date file not found: %s', xlsxfile);
    return;
end
try
    [~, ~, raw] = xlsread(xlsxfile);
catch err
    warning('Could not read %s (%s)', xlsxfile, err.message);
    return;
end
if isempty(raw), return; end

% Find header row + name / DOB columns
hdrRow = 1; nameCol = NaN; dobCol = NaN;
for r = 1:min(5, size(raw, 1))
    for c = 1:size(raw, 2)
        s = lower(strip_accents(safe_char(raw{r, c})));
        if isempty(s), continue; end
        if isnan(nameCol) && (strcmp(s, 'nom') || strcmp(s, 'name'))
            nameCol = c; hdrRow = r;
        end
        if isnan(dobCol) && (~isempty(strfind(s, 'naissance')) || ...
                ~isempty(strfind(s, 'birth')) || strcmp(s, 'dob'))
            dobCol = c; hdrRow = r;
        end
    end
    if ~isnan(nameCol) && ~isnan(dobCol), break; end
end
if isnan(nameCol) || isnan(dobCol)
    nameCol = 1; dobCol = 7; hdrRow = 1;
    warning('Birth-table headers not auto-detected; using col 1 (name) / col 7 (DOB).');
end
fprintf('Birth table: name col = %d, DOB col = %d (header row %d).\n', ...
    nameCol, dobCol, hdrRow);

% --- one-shot diagnostic: show first 3 data rows so cell types are visible
nShow = min(3, size(raw, 1) - hdrRow);
for r = hdrRow+1 : hdrRow + nShow
    if nameCol > size(raw, 2) || dobCol > size(raw, 2), continue; end
    nm = raw{r, nameCol}; bd = raw{r, dobCol};
    fprintf('  [debug] row %d: name=%s (%s) | DOB=%s (%s)\n', r, ...
        repr_cell(nm), class(nm), repr_cell(bd), class(bd));
end

nParsed = 0; nSkipName = 0; nSkipDate = 0;
for r = hdrRow+1 : size(raw, 1)
    if nameCol > size(raw, 2) || dobCol > size(raw, 2), continue; end
    nm = raw{r, nameCol};
    nmStr = to_name_string(nm);
    if isempty(nmStr), nSkipName = nSkipName + 1; continue; end
    bd = raw{r, dobCol};
    dn = parse_french_date(bd);
    if isnan(dn), nSkipDate = nSkipDate + 1; continue; end
    birthMap(lower(nmStr)) = dn;
    nParsed = nParsed + 1;
end
fprintf('Birth table: parsed %d entries (skipped %d for name, %d for DOB) from %s\n', ...
    nParsed, nSkipName, nSkipDate, xlsxfile);
end


function s = to_name_string(x)
% Normalise a cell value into a trimmed name string ('' if not usable).
if ischar(x), s = strtrim(x);
elseif isstring(x) && isscalar(x), s = strtrim(char(x));
elseif iscell(x) && isscalar(x), s = to_name_string(x{1});
elseif isnumeric(x) && isscalar(x) && ~isnan(x), s = num2str(x);
else, s = '';
end
end


function s = repr_cell(x)
% Tiny helper for the debug print.
if ischar(x), s = ['''' x ''''];
elseif isstring(x) && isscalar(x), s = ['"' char(x) '"'];
elseif isnumeric(x) && isscalar(x), s = num2str(x);
elseif isa(x,'datetime') && isscalar(x), s = datestr(x);
elseif isempty(x), s = '[]';
else, try, s = char(string(x)); catch, s = '<?>'; end
end
end


function [tdn, tval] = load_temperature_table(txtfile)
% Temperature.txt: tab-delimited
%   <idx>\t<Date YYYY-MM-DD>\t<Heure HH:MM:SS>\t<TEMP C>\t<CO>\t<alarm_hot>\t<alarm_cold>
% one header row. The file uses '---' to mark missing measurements, which
% would stop a %f format on that field and truncate the read. We therefore
% read *all* fields as strings and convert numerics with str2double, which
% turns '---' into NaN cleanly.
tdn = []; tval = [];
if exist(txtfile, 'file') ~= 2
    warning('Temperature file not found: %s', txtfile);
    return;
end
fid = fopen(txtfile, 'r');
fgetl(fid);   % skip header
C = textscan(fid, '%s%s%s%s%s%s%s', 'Delimiter', '\t');
fclose(fid);

% align column lengths (last row in some logs has fewer fields)
n = min(cellfun(@numel, C));
for k = 1:numel(C), C{k} = C{k}(1:n); end
ds   = C{2};
ts   = C{3};
tval = str2double(C{4});   % '---' -> NaN, '18.5' -> 18.5

% datetime parsing
fullStr = cell(n, 1);
for i = 1:n
    if ischar(ds{i}) && ischar(ts{i})
        fullStr{i} = [ds{i} ' ' ts{i}];
    else
        fullStr{i} = '';
    end
end
try
    tdn = datenum(fullStr, 'yyyy-mm-dd HH:MM:SS');
catch
    tdn = NaN(n, 1);
    for i = 1:n
        try, tdn(i) = datenum(fullStr{i}, 'yyyy-mm-dd HH:MM:SS'); catch, tdn(i) = NaN; end
    end
end

ok = ~isnan(tdn) & ~isnan(tval);
tdn = tdn(ok); tval = tval(ok);
[tdn, o] = sort(tdn);
tval = tval(o);

fprintf('Loaded %d temperature samples from %s\n', numel(tdn), txtfile);
if ~isempty(tdn)
    fprintf('  Temperature coverage: %s -> %s\n', ...
        datestr(tdn(1),   'yyyy-mm-dd HH:MM'), ...
        datestr(tdn(end), 'yyyy-mm-dd HH:MM'));
end
end


function dn = parse_french_date(val)
% Parse a French date. Accepts:
%   - char  : "3 juin 2020" / "03/06/2020" / "2020-06-03"
%   - string: same as char
%   - numeric scalar: Excel serial date
%   - datetime: returned directly
%   - 1-element cell: unwrap and retry
dn = NaN;
if iscell(val) && isscalar(val), val = val{1}; end
if isa(val, 'datetime') && isscalar(val) && ~isnat(val)
    dn = datenum(val); return;
end
if isnumeric(val) && isscalar(val) && ~isnan(val)
    dn = datenum('1899-12-30', 'yyyy-mm-dd') + val;   % Excel epoch
    return;
end
if isstring(val) && isscalar(val), val = char(val); end
if ~ischar(val), return; end
s = strtrim(val);
if isempty(s), return; end

% (1) French month name: "3 juin 2020"
s2 = strip_accents(lower(s));
tok = regexp(s2, '(\d{1,2})\s+([a-z]+)\.?\s+(\d{4})', 'tokens', 'once');
if ~isempty(tok)
    months = {'janvier','fevrier','mars','avril','mai','juin', ...
              'juillet','aout','septembre','octobre','novembre','decembre'};
    % also handle 3-letter abbreviations (janv., fevr., juil., aout, sept., oct., nov., dec.)
    abbr3  = {'jan','fev','mar','avr','mai','jui','jui','aou','sep','oct','nov','dec'};
    day_ = str2double(tok{1}); yr_ = str2double(tok{3}); mon_ = tok{2};
    midx = find(strcmp(mon_, months), 1);
    if isempty(midx)
        midx = find(strncmp(mon_, abbr3, 3), 1);
        % distinguish 'jui' (June vs July): require >=4 chars to disambiguate
        if ~isempty(midx) && strncmp(mon_,'jui',3)
            if length(mon_) >= 4 && mon_(4) == 'l'
                midx = 7;   % juillet
            else
                midx = 6;   % juin
            end
        end
    end
    if ~isempty(midx)
        dn = datenum(yr_, midx, day_); return;
    end
end

% (2) Numeric formats: try a few common ones
fmts = {'dd/mm/yyyy','d/m/yyyy','dd-mm-yyyy','yyyy-mm-dd','yyyy/mm/dd','dd.mm.yyyy'};
for k = 1:numel(fmts)
    try, dn = datenum(s, fmts{k}); return; catch, end %#ok<NOSEMI>
end
end


function s = strip_accents(s)
if ~ischar(s), s = ''; return; end
from = {char(233), char(232), char(234), char(235), ...   % e
        char(224), char(226), ...                         % a
        char(244), char(238), char(239), char(251), char(249), char(231)};
to   = {'e','e','e','e','a','a','o','i','i','u','u','c'};
for k = 1:numel(from)
    s = strrep(s, from{k}, to{k});
end
end


function s = safe_char(x)
if ischar(x), s = x;
elseif isnumeric(x) && isscalar(x), s = num2str(x);
else, s = '';
end
end


function out = sec2hms(sec)
if isnan(sec), out = ''; return; end
sec = round(sec);
h = floor(sec / 3600);
m = floor(mod(sec, 3600) / 60);
s = mod(sec, 60);
out = sprintf('%02d:%02d:%02d', h, m, s);
end


function out = datestr_safe(dt)
if isa(dt, 'datetime') && ~isnat(dt)
    out = datestr(dt, 'HH:MM:SS');
else
    out = '--:--:--';
end
end


function out = round_or_nan(x)
if isnan(x), out = -1; else, out = round(x); end
end


function hms = extract_hms(v)
% Extract [H M S] from [H M S] or [Y M D H M S].
hms = [NaN NaN NaN];
v = v(:).';
if numel(v) == 3
    hms = v;
elseif numel(v) >= 6
    hms = v(4:6);
end
end