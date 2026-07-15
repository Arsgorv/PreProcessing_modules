function session = make_ns2_metadata_from_oe(datapath)

if nargin < 1 || isempty(datapath)
    datapath = uigetdir(pwd,'Select Open Ephys continuous stream folder');
end

datapath = char(datapath);
datFile = fullfile(datapath,'continuous.dat');

if ~exist(datFile,'file')
    error('continuous.dat not found in: %s',datapath);
end

basename = 'continuous';

structureFile = find_file_upwards(datapath,'structure.oebin',12);
if isempty(structureFile)
    error('structure.oebin not found above datapath');
end

settingsFile = find_settings_xml(datapath,12);

fprintf('structure.oebin: %s\n',structureFile);
if isempty(settingsFile)
    warning('settings.xml not found. Continuing without it.');
else
    fprintf('settings.xml: %s\n',settingsFile);
end

txt = fileread(structureFile);
oe = jsondecode(txt);

if ~isfield(oe,'continuous')
    error('No continuous field found in structure.oebin');
end

continuousList = oe.continuous;

recordingPath = fileparts(structureFile);
relativeFolder = strrep(datapath,[recordingPath filesep],'');
relativeFolder = normalize_oe_path(relativeFolder);

leafFolder = get_last_folder_name(datapath);
leafFolderNorm = normalize_oe_path(leafFolder);

nStreams = length(continuousList);
scores = zeros(1,nStreams);
streamInfo = cell(nStreams,1);

fprintf('\nContinuous streams in structure.oebin:\n');

for i = 1:nStreams
    entry = get_entry(continuousList,i);

    folderName = get_text_field(entry,'folder_name');
    streamName = get_text_field(entry,'stream_name');
    sourceName = get_text_field(entry,'source_processor_name');

    folderNorm = normalize_oe_path(folderName);
    streamNorm = normalize_oe_path(streamName);
    sourceNorm = normalize_oe_path(sourceName);

    nCh = get_num_field(entry,'num_channels');
    sr = get_num_field(entry,'sample_rate');

    score = 0;

    if strcmp(remove_slashes(folderNorm),remove_slashes(leafFolderNorm))
        score = score + 200;
    end

    if strcmp(remove_slashes(streamNorm),remove_slashes(leafFolderNorm))
        score = score + 150;
    end

    if strcmp(folderNorm,relativeFolder)
        score = score + 100;
    end

    if has_substring(relativeFolder,folderNorm) || has_substring(folderNorm,relativeFolder)
        score = score + 50;
    end

    if has_substring(folderNorm,leafFolderNorm) || has_substring(leafFolderNorm,folderNorm)
        score = score + 40;
    end

    if has_substring(streamNorm,leafFolderNorm) || has_substring(leafFolderNorm,streamNorm)
        score = score + 30;
    end

    if has_substring(sourceNorm,leafFolderNorm) || has_substring(leafFolderNorm,sourceNorm)
        score = score + 20;
    end

    if has_substring(leafFolderNorm,'lfp') && has_substring(folderNorm,'lfp')
        score = score + 100;
    end

    if has_substring(leafFolderNorm,'lfp') && has_substring(streamNorm,'lfp')
        score = score + 100;
    end

    if has_substring(leafFolderNorm,'lfp') && sr > 2000 && sr < 4000
        score = score + 60;
    end

    if has_substring(leafFolderNorm,'ap') && has_substring(folderNorm,'ap')
        score = score + 100;
    end

    if has_substring(leafFolderNorm,'ap') && has_substring(streamNorm,'ap')
        score = score + 100;
    end

    if has_substring(leafFolderNorm,'ap') && sr > 20000
        score = score + 60;
    end

    if has_substring(leafFolderNorm,'adc') && has_substring(folderNorm,'adc')
        score = score + 100;
    end

    scores(i) = score;
    streamInfo{i} = folderName;

    fprintf('%2d | score %3d | sr %8.1f | nCh %4.0f | folder: %s | stream: %s\n', ...
        i, score, sr, nCh, folderName, streamName);
end

[bestScore,bestIdx] = max(scores);

if bestScore == 0
    error('Could not match datapath to a stream in structure.oebin. Check printed list above.');
end

entry = get_entry(continuousList,bestIdx);

nChannels = get_num_field(entry,'num_channels');
sampleRate = get_num_field(entry,'sample_rate');

if isnan(nChannels) || nChannels <= 0
    if isfield(entry,'channels')
        nChannels = length(entry.channels);
    end
end

if isnan(sampleRate) || sampleRate <= 0
    if has_substring(leafFolderNorm,'lfp')
        sampleRate = 2500;
        warning('sample_rate missing. Forced sampleRate = 2500 for LFP.');
    elseif has_substring(leafFolderNorm,'ap')
        sampleRate = 30000;
        warning('sample_rate missing. Forced sampleRate = 30000 for AP.');
    else
        error('sample_rate missing and stream type unclear.');
    end
end

if isnan(nChannels) || nChannels <= 0
    if has_substring(leafFolderNorm,'lfp') || has_substring(leafFolderNorm,'ap')
        nChannels = 384;
        warning('num_channels missing. Forced nChannels = 384.');
    else
        error('num_channels missing and stream type unclear.');
    end
end

lsb = get_bit_volts(entry);
if isnan(lsb)
    lsb = 0.195;
    warning('bit_volts missing. Forced leastSignificantBit = 0.195.');
end

fprintf('\nSelected stream %d:\n',bestIdx);
fprintf('folder = %s\n',streamInfo{bestIdx});
fprintf('nChannels = %d\n',nChannels);
fprintf('sampleRate = %.1f Hz\n',sampleRate);
fprintf('leastSignificantBit = %.6f\n',lsb);

session = struct;

session.general.name = basename;
session.general.basePath = datapath;

session.extracellular.fileName = 'continuous.dat';
session.extracellular.nChannels = nChannels;
session.extracellular.sr = sampleRate;
session.extracellular.srLfp = sampleRate;
session.extracellular.precision = 'int16';
session.extracellular.leastSignificantBit = lsb;

session.extracellular.nElectrodeGroups = 1;
session.extracellular.electrodeGroups.channels = {1:nChannels};

session.general.openEphys.structureFile = structureFile;
session.general.openEphys.settingsFile = settingsFile;
session.general.openEphys.selectedStream = bestIdx;
session.general.openEphys.selectedStreamFolder = streamInfo{bestIdx};

outFile = fullfile(datapath,[basename '.session.mat']);
save(outFile,'session');

fprintf('\nSaved: %s\n',outFile);

end


function filePath = find_file_upwards(startPath,fileName,maxDepth)

filePath = '';
searchPath = startPath;

for k = 1:maxDepth
    candidate = fullfile(searchPath,fileName);

    if exist(candidate,'file')
        filePath = candidate;
        return
    end

    parentPath = fileparts(searchPath);

    if strcmp(parentPath,searchPath)
        return
    end

    searchPath = parentPath;
end

end


function settingsFile = find_settings_xml(startPath,maxDepth)

settingsFile = find_file_upwards(startPath,'settings.xml',maxDepth);

if ~isempty(settingsFile)
    return
end

searchPath = startPath;

for k = 1:maxDepth
    d = dir(fullfile(searchPath,'Record Node*'));

    for i = 1:length(d)
        if d(i).isdir
            candidate = fullfile(searchPath,d(i).name,'settings.xml');

            if exist(candidate,'file')
                settingsFile = candidate;
                return
            end
        end
    end

    parentPath = fileparts(searchPath);

    if strcmp(parentPath,searchPath)
        return
    end

    searchPath = parentPath;
end

end


function entry = get_entry(list,i)

if iscell(list)
    entry = list{i};
else
    entry = list(i);
end

end


function txt = get_text_field(s,fieldName)

txt = '';

if isfield(s,fieldName)
    val = s.(fieldName);

    if ischar(val)
        txt = val;
    elseif isstring(val)
        txt = char(val);
    elseif isnumeric(val)
        txt = num2str(val);
    end
end

end


function x = get_num_field(s,fieldName)

x = NaN;

if isfield(s,fieldName)
    val = s.(fieldName);

    if isnumeric(val)
        x = double(val);
    elseif ischar(val)
        x = str2double(val);
    elseif isstring(val)
        x = str2double(char(val));
    end
end

end


function p = normalize_oe_path(p)

p = char(p);
p = strrep(p,'\','/');
p = strrep(p,'//','/');
p = lower(p);

while ~isempty(p) && p(end) == '/'
    p(end) = [];
end

end


function p = remove_slashes(p)

p = strrep(p,'/','');
p = strrep(p,'\','');

end


function tf = has_substring(a,b)

tf = false;

if isempty(a) || isempty(b)
    return
end

tf = ~isempty(strfind(a,b));

end


function folderName = get_last_folder_name(pathName)

pathName = char(pathName);

while ~isempty(pathName) && (pathName(end) == filesep || pathName(end) == '/' || pathName(end) == '\')
    pathName(end) = [];
end

pathName = strrep(pathName,'\','/');
parts = strsplit(pathName,'/');
folderName = parts{end};

end


function lsb = get_bit_volts(entry)

lsb = NaN;

if ~isfield(entry,'channels')
    return
end

try
    ch = entry.channels;
    bitvolts = [ch.bit_volts];
    bitvolts = double(bitvolts);
    bitvolts = bitvolts(~isnan(bitvolts));

    if ~isempty(bitvolts)
        lsb = median(bitvolts);
    end
catch
    lsb = NaN;
end

end