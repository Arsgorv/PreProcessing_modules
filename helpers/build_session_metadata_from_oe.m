function session_metadata = build_session_metadata_from_oe(sessions)
% build_session_metadata_from_oe_AG
%
% Build per-session metadata from Open Ephys / ExpeInfo files.
%
% Input:
%   sessions
%       cell array, string array, char, table, or struct containing session paths.
%
% Output:
%   session_metadata
%       table with recording date/start/end, duration, animal, and placeholders
%       for animal weight and age.
%
% Main logic:
%   1. Find Open Ephys recording folders.
%   2. Get absolute start time, preferably from ExpeInfo.mat.
%   3. If unavailable, try sync_messages.txt/settings.xml/folder names.
%   4. Get duration from continuous.dat size and structure.oebin.
%   5. Compute recording_end = recording_start + duration when possible.
%
% Usage:
%   session_metadata = build_session_metadata_from_oe_AG(sessions);

nSessions = get_num_sessions_AG(sessions);

Session = cell(nSessions,1);
Animal = cell(nSessions,1);
Datapath = cell(nSessions,1);

RecordingStart = NaT(nSessions,1);
RecordingEnd = NaT(nSessions,1);
RecordingDate = NaT(nSessions,1);

RecordingDateSource = cell(nSessions,1);
RecordingStartSource = cell(nSessions,1);

SessionDurationSec = nan(nSessions,1);
TotalRecordedDurationSec = nan(nSessions,1);
NOERecordings = nan(nSessions,1);

HasExpeInfo = false(nSessions,1);
HasSettingsXML = false(nSessions,1);
HasSyncMessages = false(nSessions,1);

AnimalWeight_g = nan(nSessions,1);
AnimalAge_days = nan(nSessions,1);

OERecordingFolders = cell(nSessions,1);
ParserNotes = cell(nSessions,1);

for iSess = 1:nSessions

    datapath = get_session_path_AG(sessions, iSess);
    datapath = remove_trailing_filesep_AG(datapath);

    Datapath{iSess} = datapath;

    [~, sessionName] = fileparts(datapath);
    Session{iSess} = sessionName;

    Animal{iSess} = derive_animal_name_from_path_AG(datapath);

    HasExpeInfo(iSess) = exist(fullfile(datapath, 'ExpeInfo.mat'), 'file') == 2;
    HasSettingsXML(iSess) = ~isempty(find_files_recursive_AG(datapath, 'settings.xml'));
    HasSyncMessages(iSess) = ~isempty(find_files_recursive_AG(datapath, 'sync_messages.txt'));

    [recStart, recEnd, recDurSec, oeFolders, notes, startSource] = summarize_oe_recordings_AG(datapath);

    NOERecordings(iSess) = numel(oeFolders);
    OERecordingFolders{iSess} = strjoin(oeFolders, ' | ');
    ParserNotes{iSess} = strjoin(notes, ' | ');

    if any(~isnat(recStart))
        RecordingStart(iSess) = min(recStart(~isnat(recStart)));
        RecordingStartSource{iSess} = startSource;

        if any(~isnat(recEnd))
            RecordingEnd(iSess) = max(recEnd(~isnat(recEnd)));
        end

        RecordingDate(iSess) = dateshift(RecordingStart(iSess), 'start', 'day');
        RecordingDateSource{iSess} = 'derived from recording_start';

    else
        [dateFromPath, dateNote] = parse_date_from_path_AG(datapath);

        if ~isnat(dateFromPath)
            RecordingDate(iSess) = dateFromPath;
            RecordingDateSource{iSess} = ['session path: ' dateNote];
        else
            RecordingDateSource{iSess} = 'date not found';
        end

        RecordingStartSource{iSess} = 'exact clock time not found';
    end

    if ~isnat(RecordingStart(iSess)) && ~isnat(RecordingEnd(iSess))
        SessionDurationSec(iSess) = seconds(RecordingEnd(iSess) - RecordingStart(iSess));
    end

    if any(~isnan(recDurSec))
        TotalRecordedDurationSec(iSess) = sum(recDurSec(~isnan(recDurSec)));

        if ~isnat(RecordingStart(iSess)) && isnat(RecordingEnd(iSess))
            RecordingEnd(iSess) = RecordingStart(iSess) + seconds(TotalRecordedDurationSec(iSess));
            SessionDurationSec(iSess) = TotalRecordedDurationSec(iSess);
        end
    end

end

RecordingStart.Format = 'yyyy-MM-dd HH:mm:ss';
RecordingEnd.Format = 'yyyy-MM-dd HH:mm:ss';
RecordingDate.Format = 'yyyy-MM-dd';

session_metadata = table( ...
    Session, ...
    Animal, ...
    Datapath, ...
    RecordingStart, ...
    RecordingEnd, ...
    RecordingDate, ...
    RecordingDateSource, ...
    RecordingStartSource, ...
    SessionDurationSec, ...
    TotalRecordedDurationSec, ...
    NOERecordings, ...
    HasExpeInfo, ...
    HasSettingsXML, ...
    HasSyncMessages, ...
    AnimalWeight_g, ...
    AnimalAge_days, ...
    OERecordingFolders, ...
    ParserNotes);

session_metadata.Properties.VariableNames = { ...
    'session', ...
    'animal', ...
    'datapath', ...
    'recording_start', ...
    'recording_end', ...
    'recording_date', ...
    'recording_date_source', ...
    'recording_start_source', ...
    'session_duration_sec', ...
    'total_recorded_duration_sec', ...
    'n_oe_recordings', ...
    'has_expeinfo', ...
    'has_settings_xml', ...
    'has_sync_messages', ...
    'animal_weight_g', ...
    'animal_age_days', ...
    'oe_recording_folders', ...
    'parser_notes'};

end


function [recStart, recEnd, recDurSec, oeFolders, notes, startSource] = summarize_oe_recordings_AG(datapath)

oeFolders = find_oe_recording_folders_AG(datapath);

nRec = numel(oeFolders);

recStart = NaT(nRec,1);
recEnd = NaT(nRec,1);
recDurSec = nan(nRec,1);
notes = cell(nRec,1);
startSource = '';

[expeStart, expeNote] = parse_start_datetime_from_expeinfo_AG(datapath);

if nRec == 0
    notes = {'No Open Ephys recording folders found'};
    startSource = expeNote;
    return
end

for iRec = 1:nRec

    oeFolder = oeFolders{iRec};

    [thisStart, startNote] = parse_oe_start_datetime_AG(oeFolder);

    if isnat(thisStart) && ~isnat(expeStart)
        thisStart = expeStart;
        startNote = expeNote;
    end

    [thisDurSec, durNote] = parse_oe_duration_seconds_AG(oeFolder);

    recStart(iRec) = thisStart;
    recDurSec(iRec) = thisDurSec;

    if ~isnat(thisStart) && ~isnan(thisDurSec)
        recEnd(iRec) = thisStart + seconds(thisDurSec);
    end

    notes{iRec} = [startNote ', ' durNote];

end

validStarts = recStart(~isnat(recStart));

if ~isempty(validStarts)
    startSource = notes{find(~isnat(recStart), 1)};
else
    startSource = 'exact clock time not found';
end

end


function oeFolders = find_oe_recording_folders_AG(datapath)

oeFolders = {};

oebinFiles = find_files_recursive_AG(datapath, 'structure.oebin');
syncFiles = find_files_recursive_AG(datapath, 'sync_messages.txt');
oldContinuousFiles = find_files_recursive_AG(datapath, '*.continuous');

for iFile = 1:numel(oebinFiles)
    oeFolders{end+1,1} = fileparts(oebinFiles{iFile});
end

for iFile = 1:numel(syncFiles)
    oeFolders{end+1,1} = fileparts(syncFiles{iFile});
end

for iFile = 1:numel(oldContinuousFiles)
    oeFolders{end+1,1} = fileparts(oldContinuousFiles{iFile});
end

if isempty(oeFolders)
    return
end

oeFolders = unique(oeFolders, 'stable');

end


function [tStart, note] = parse_start_datetime_from_expeinfo_AG(datapath)

tStart = NaT;
note = 'start datetime not found in ExpeInfo';

expeFile = fullfile(datapath, 'ExpeInfo.mat');

if ~exist(expeFile, 'file')
    note = 'ExpeInfo.mat not found';
    return
end

try
    S = load(expeFile);
catch
    note = 'ExpeInfo.mat found but could not be loaded';
    return
end

if ~isfield(S, 'ExpeInfo')
    note = 'ExpeInfo.mat found but no ExpeInfo variable';
    return
end

ExpeInfo = S.ExpeInfo;

candidateText = {};

if isfield(ExpeInfo, 'PreProcessingInfo')

    P = ExpeInfo.PreProcessingInfo;

    if isfield(P, 'FolderForConcatenation_Ephys')
        candidateText = add_text_candidate_AG(candidateText, P.FolderForConcatenation_Ephys);
    end

    if isfield(P, 'FolderSessionName')
        candidateText = add_text_candidate_AG(candidateText, P.FolderSessionName);
    end

end

if isfield(ExpeInfo, 'SessionType')
    candidateText = add_text_candidate_AG(candidateText, ExpeInfo.SessionType);
end

if isfield(ExpeInfo, 'date')
    candidateText = add_text_candidate_AG(candidateText, ExpeInfo.date);
end

for iText = 1:numel(candidateText)

    thisText = candidateText{iText};

    [thisTime, thisNote] = parse_datetime_from_text_AG(thisText);

    if ~isnat(thisTime)
        tStart = thisTime;
        note = ['start datetime from ExpeInfo: ' thisNote];
        return
    end

end

for iText = 1:numel(candidateText)

    thisText = candidateText{iText};

    [thisDate, thisNote] = parse_date_from_text_AG(thisText);

    if ~isnat(thisDate)
        tStart = thisDate;
        note = ['date only from ExpeInfo: ' thisNote];
        return
    end

end

end


function candidateText = add_text_candidate_AG(candidateText, x)

if ischar(x)
    candidateText{end+1,1} = x;

elseif isstring(x)
    candidateText{end+1,1} = char(x);

elseif iscell(x)
    for i = 1:numel(x)
        candidateText = add_text_candidate_AG(candidateText, x{i});
    end

end

end


function [tStart, note] = parse_oe_start_datetime_AG(oeFolder)

tStart = NaT;
note = 'start datetime not found';

syncFile = fullfile(oeFolder, 'sync_messages.txt');

if exist(syncFile, 'file')
    txt = fileread(syncFile);

    [tStart, okNote] = parse_datetime_from_text_AG(txt);

    if ~isnat(tStart)
        note = ['start from sync_messages text: ' okNote];
        return
    end

    [tStart, okNote] = parse_software_time_AG(txt);

    if ~isnat(tStart)
        note = ['start from sync_messages software time: ' okNote];
        return
    end
end

settingsFile = fullfile(oeFolder, 'settings.xml');

if exist(settingsFile, 'file')
    txt = fileread(settingsFile);

    [tStart, okNote] = parse_datetime_from_text_AG(txt);

    if ~isnat(tStart)
        note = ['start from settings.xml: ' okNote];
        return
    end
end

[tStart, okNote] = parse_datetime_from_path_AG(oeFolder);

if ~isnat(tStart)
    note = ['start from folder name: ' okNote];
    return
end

[tDate, okNote] = parse_date_from_path_AG(oeFolder);

if ~isnat(tDate)
    tStart = tDate;
    note = ['date only from folder name: ' okNote];
end

end


function [durSec, note] = parse_oe_duration_seconds_AG(oeFolder)

durSec = nan;
note = 'duration not found';

oebinFile = fullfile(oeFolder, 'structure.oebin');

if exist(oebinFile, 'file')
    [durSec, okNote] = duration_from_oebin_AG(oeFolder, oebinFile);

    if ~isnan(durSec)
        note = ['duration from structure.oebin/continuous.dat: ' okNote];
        return
    end
end

[durSec, okNote] = duration_from_old_continuous_AG(oeFolder);

if ~isnan(durSec)
    note = ['duration from old .continuous file: ' okNote];
end

end


function [durSec, note] = duration_from_oebin_AG(oeFolder, oebinFile)

durSec = nan;
note = 'failed';

try
    info = jsondecode(fileread(oebinFile));
catch
    note = 'could not jsondecode structure.oebin';
    return
end

if ~isfield(info, 'continuous')
    note = 'no continuous field in structure.oebin';
    return
end

continuousInfo = info.continuous;

durList = [];

for iCont = 1:numel(continuousInfo)

    thisInfo = continuousInfo(iCont);

    fs = nan;
    nChan = nan;
    folderName = '';

    if isfield(thisInfo, 'sample_rate')
        fs = double(thisInfo.sample_rate);
    elseif isfield(thisInfo, 'sampleRate')
        fs = double(thisInfo.sampleRate);
    end

    if isfield(thisInfo, 'num_channels')
        nChan = double(thisInfo.num_channels);
    elseif isfield(thisInfo, 'numChannels')
        nChan = double(thisInfo.numChannels);
    elseif isfield(thisInfo, 'channels')
        nChan = numel(thisInfo.channels);
    end

    if isfield(thisInfo, 'folder_name')
        folderName = thisInfo.folder_name;
    elseif isfield(thisInfo, 'folderName')
        folderName = thisInfo.folderName;
    end

    if isnan(fs) || isnan(nChan) || nChan == 0
        continue
    end

    candidateFiles = {};

    if ~isempty(folderName)
        candidateFiles{end+1,1} = fullfile(oeFolder, 'continuous', folderName, 'continuous.dat');
        candidateFiles{end+1,1} = fullfile(oeFolder, folderName, 'continuous.dat');
    end

    datFiles = find_files_recursive_AG(oeFolder, 'continuous.dat');

    for iFile = 1:numel(datFiles)
        candidateFiles{end+1,1} = datFiles{iFile};
    end

    candidateFiles = unique(candidateFiles, 'stable');

    for iFile = 1:numel(candidateFiles)

        datFile = candidateFiles{iFile};

        if ~exist(datFile, 'file')
            continue
        end

        d = dir(datFile);

        nSamples = double(d.bytes) / 2 / nChan;   % int16 = 2 bytes/sample
        thisDurSec = nSamples / fs;

        if thisDurSec > 0
            durList(end+1,1) = thisDurSec;
        end

    end

end

if ~isempty(durList)
    durSec = max(durList);
    note = sprintf('max continuous stream %.1f sec', durSec);
end

end


function [durSec, note] = duration_from_old_continuous_AG(oeFolder)

durSec = nan;
note = 'no old .continuous files';

contFiles = find_files_recursive_AG(oeFolder, '*.continuous');

if isempty(contFiles)
    return
end

durList = nan(numel(contFiles),1);

for iFile = 1:numel(contFiles)

    thisFile = contFiles{iFile};

    fid = fopen(thisFile, 'r');

    if fid < 0
        continue
    end

    header = fread(fid, 1024, '*char')';
    fclose(fid);

    fs = parse_sample_rate_from_header_AG(header);

    if isnan(fs)
        continue
    end

    d = dir(thisFile);

    bytesPerRecord = 8 + 2 + 2 + 1024*2 + 10;
    nRecords = floor((double(d.bytes) - 1024) / bytesPerRecord);

    if nRecords > 0
        durList(iFile) = nRecords * 1024 / fs;
    end

end

if any(~isnan(durList))
    durSec = max(durList(~isnan(durList)));
    note = sprintf('max old continuous stream %.1f sec', durSec);
end

end


function fs = parse_sample_rate_from_header_AG(header)

fs = nan;

tok = regexp(header, 'sampleRate\s*=\s*([0-9\.]+)', 'tokens', 'once');

if isempty(tok)
    tok = regexp(header, 'sample_rate\s*=\s*([0-9\.]+)', 'tokens', 'once');
end

if ~isempty(tok)
    fs = str2double(tok{1});
end

end


function [t, note] = parse_datetime_from_text_AG(txt)

t = NaT;
note = '';

if ~ischar(txt)
    return
end

pattern = '(\d{4})[-_](\d{2})[-_](\d{2})[ T_](\d{2})[:\-](\d{2})[:\-](\d{2})(?:\.(\d+))?';

tok = regexp(txt, pattern, 'tokens', 'once');

if isempty(tok)
    return
end

yearVal = str2double(tok{1});
monthVal = str2double(tok{2});
dayVal = str2double(tok{3});
hourVal = str2double(tok{4});
minVal = str2double(tok{5});
secVal = str2double(tok{6});

if numel(tok) >= 7 && ~isempty(tok{7})
    fracVal = str2double(['0.' tok{7}]);
else
    fracVal = 0;
end

if isnan(yearVal) || isnan(monthVal) || isnan(dayVal) || isnan(hourVal) || isnan(minVal) || isnan(secVal)
    return
end

t = datetime(yearVal, monthVal, dayVal, hourVal, minVal, secVal + fracVal);
note = datestr(t, 'yyyy-mm-dd HH:MM:SS');

end


function [t, note] = parse_date_from_text_AG(txt)

t = NaT;
note = '';

if ~ischar(txt)
    return
end

tok = regexp(txt, '(\d{4})[-_]?(\d{2})[-_]?(\d{2})', 'tokens', 'once');

if isempty(tok)
    return
end

y = str2double(tok{1});
m = str2double(tok{2});
d = str2double(tok{3});

if isnan(y) || isnan(m) || isnan(d)
    return
end

if y < 2010 || y > 2035 || m < 1 || m > 12 || d < 1 || d > 31
    return
end

t = datetime(y, m, d);
note = datestr(t, 'yyyy-mm-dd');

end


function [t, note] = parse_software_time_AG(txt)

t = NaT;
note = '';

tok = regexp(txt, 'Software time:\s*([0-9]+)', 'tokens', 'once');

if isempty(tok)
    return
end

x = str2double(tok{1});

if isnan(x)
    return
end

if x > 1e15
    posixTime = x / 1e6;
    unitNote = 'microseconds';
elseif x > 1e12
    posixTime = x / 1e3;
    unitNote = 'milliseconds';
elseif x > 1e9
    posixTime = x;
    unitNote = 'seconds';
else
    return
end

try
    tCandidate = datetime(posixTime, 'ConvertFrom', 'posixtime', 'TimeZone', 'local');
    tCandidate.TimeZone = '';

    yr = year(tCandidate);

    if yr >= 2010 && yr <= 2035
        t = tCandidate;
        note = ['Software time interpreted as Unix ' unitNote];
    end
catch
    t = NaT;
end

end


function [t, note] = parse_datetime_from_path_AG(pathStr)

t = NaT;
note = '';

parts = regexp(pathStr, '[\\/]', 'split');
parts = parts(~cellfun('isempty', parts));

for iPart = numel(parts):-1:1

    thisPart = parts{iPart};

    [thisTime, thisNote] = parse_datetime_from_text_AG(thisPart);

    if ~isnat(thisTime)
        t = thisTime;
        note = ['path component ' thisPart ', ' thisNote];
        return
    end

end

end


function [dateValue, note] = parse_date_from_path_AG(pathStr)

dateValue = NaT;
note = '';

parts = regexp(pathStr, '[\\/]', 'split');
parts = parts(~cellfun('isempty', parts));

for iPart = numel(parts):-1:1

    thisPart = parts{iPart};

    [thisDate, thisNote] = parse_date_from_text_AG(thisPart);

    if ~isnat(thisDate)
        dateValue = thisDate;
        note = ['path component ' thisPart ', ' thisNote];
        return
    end

end

end


function animal = derive_animal_name_from_path_AG(datapath)

datapath = remove_trailing_filesep_AG(datapath);

parts = regexp(datapath, '[\\/]', 'split');
parts = parts(~cellfun('isempty', parts));

animal = '';

if isempty(parts)
    return
end

dateIdx = [];

for iPart = 1:numel(parts)
    if is_date_like_component_AG(parts{iPart})
        dateIdx = iPart;
    end
end

phaseNames = { ...
    'freely-moving', ...
    'free-moving', ...
    'head-fixed', ...
    'headfixed', ...
    'hf', ...
    'fm', ...
    'PreSleep', ...
    'PostSleep', ...
    'Conditioning', ...
    'PostTest', ...
    'sleep', ...
    'wake'};

if ~isempty(dateIdx)

    if dateIdx > 2 && any(strcmpi(parts{dateIdx - 1}, phaseNames))
        animal = parts{dateIdx - 2};
        return
    end

    if dateIdx > 1
        animal = parts{dateIdx - 1};
        return
    end

end

if numel(parts) >= 2
    animal = parts{end - 1};
else
    animal = parts{end};
end

end


function tf = is_date_like_component_AG(s)

tf = false;

if ~isempty(regexp(s, '^\d{4}[-_]\d{2}[-_]\d{2}', 'once'))
    tf = true;
    return
end

if ~isempty(regexp(s, '^\d{8}', 'once'))
    tf = true;
    return
end

end


function files = find_files_recursive_AG(rootFolder, pattern)

files = {};

if ~exist(rootFolder, 'dir')
    return
end

dFiles = dir(fullfile(rootFolder, pattern));

for iFile = 1:numel(dFiles)
    if ~dFiles(iFile).isdir
        files{end+1,1} = fullfile(rootFolder, dFiles(iFile).name);
    end
end

d = dir(rootFolder);

for i = 1:numel(d)

    if ~d(i).isdir
        continue
    end

    if strcmp(d(i).name, '.') || strcmp(d(i).name, '..')
        continue
    end

    subFolder = fullfile(rootFolder, d(i).name);
    subFiles = find_files_recursive_AG(subFolder, pattern);

    for j = 1:numel(subFiles)
        files{end+1,1} = subFiles{j};
    end

end

end


function nSessions = get_num_sessions_AG(sessions)

if ischar(sessions)
    nSessions = 1;

elseif isstring(sessions)
    nSessions = numel(sessions);

elseif iscell(sessions)
    nSessions = numel(sessions);

elseif istable(sessions)
    nSessions = height(sessions);

elseif isstruct(sessions)
    nSessions = numel(sessions);

else
    error('Unsupported sessions variable type.');
end

end


function datapath = get_session_path_AG(sessions, iSess)

if ischar(sessions)
    datapath = sessions;
    return
end

if isstring(sessions)
    datapath = char(sessions(iSess));
    return
end

if iscell(sessions)

    datapath = sessions{iSess};

    if isstring(datapath)
        datapath = char(datapath);
    end

    return
end

candidateFields = {'datapath', 'DataPath', 'path', 'Path', 'folder', 'Folder', 'SessionPath'};

if istable(sessions)

    varNames = sessions.Properties.VariableNames;

    idx = [];

    for iField = 1:numel(candidateFields)

        thisIdx = find(strcmp(varNames, candidateFields{iField}), 1);

        if ~isempty(thisIdx)
            idx = thisIdx;
            break
        end

    end

    if isempty(idx)
        error('sessions table must contain a datapath/path/folder variable.');
    end

    datapath = sessions{iSess, idx};

    if iscell(datapath)
        datapath = datapath{1};
    end

    if isstring(datapath)
        datapath = char(datapath);
    end

    return
end

if isstruct(sessions)

    fieldName = '';

    for iField = 1:numel(candidateFields)

        if isfield(sessions, candidateFields{iField})
            fieldName = candidateFields{iField};
            break
        end

    end

    if isempty(fieldName)
        error('sessions struct must contain a datapath/path/folder field.');
    end

    datapath = sessions(iSess).(fieldName);

    if isstring(datapath)
        datapath = char(datapath);
    end

    return
end

error('Unsupported sessions variable type.');

end


function p = remove_trailing_filesep_AG(p)

if isstring(p)
    p = char(p);
end

while ~isempty(p) && (strcmp(p(end), filesep) || strcmp(p(end), '/') || strcmp(p(end), '\'))
    p(end) = [];
end

end