%% Copy all Open Ephys settings.xml files from one animal folder
% MATLAB R2018b compatible

clear
clc

sourceRoot = 'F:\OB_Sleep_project\Labneh\head-fixed';
destRoot   = 'E:\Labneh\head-fixed';

dryRun = false;   % true = only print what would be copied, false = actually copy

if ~exist(sourceRoot, 'dir')
    error('Source folder does not exist: %s', sourceRoot)
end

if ~exist(destRoot, 'dir')
    mkdir(destRoot)
end

settingsFiles = dir(fullfile(sourceRoot, '**', 'settings.xml'));

fprintf('Found %d settings.xml files\n\n', numel(settingsFiles))

copyLog = cell(numel(settingsFiles), 5);
nCopied = 0;
nSkipped = 0;

for i = 1:numel(settingsFiles)

    srcFile = fullfile(settingsFiles(i).folder, settingsFiles(i).name);

    % Usually:
    % ...\recording_folder\Record Node 101\settings.xml
    % So the relevant name is the folder above "Record Node 101"
    [settingsParent, ~, ~] = fileparts(srcFile);
    [parentOfSettingsParent, settingsParentName, ~] = fileparts(settingsParent);

    if contains(lower(settingsParentName), 'record node')
        [~, recordingFolderName, ~] = fileparts(parentOfSettingsParent);
    else
        recordingFolderName = settingsParentName;
    end

    recordingFolderName = sanitize_folder_name(recordingFolderName);

    destFolder = fullfile(destRoot, recordingFolderName);

    % Avoid overwriting if the same recording folder name appears twice
    if exist(fullfile(destFolder, 'settings.xml'), 'file')
        sessionName = get_session_folder_name(srcFile, sourceRoot);
        destFolder = fullfile(destRoot, [sessionName '__' recordingFolderName]);
        destFolder = make_unique_folder_if_needed(destFolder);
    end

    destFile = fullfile(destFolder, 'settings.xml');

    copyLog{i, 1} = srcFile;
    copyLog{i, 2} = destFile;
    copyLog{i, 3} = recordingFolderName;
    copyLog{i, 4} = settingsFiles(i).bytes;
    copyLog{i, 5} = '';

    fprintf('[%d/%d]\n', i, numel(settingsFiles))
    fprintf('FROM: %s\n', srcFile)
    fprintf('TO:   %s\n', destFile)

    if dryRun
        fprintf('DRY RUN: not copied\n\n')
        copyLog{i, 5} = 'dry run';
        continue
    end

    if ~exist(destFolder, 'dir')
        mkdir(destFolder)
    end

    [ok, msg] = copyfile(srcFile, destFile);

    if ok
        fprintf('Copied\n\n')
        copyLog{i, 5} = 'copied';
        nCopied = nCopied + 1;
    else
        fprintf('FAILED: %s\n\n', msg)
        copyLog{i, 5} = ['failed: ' msg];
        nSkipped = nSkipped + 1;
    end
end

logTable = cell2table(copyLog, ...
    'VariableNames', {'source_file', 'destination_file', 'recording_folder', 'size_bytes', 'status'});

logFile = fullfile(destRoot, 'settings_xml_copy_log.csv');
writetable(logTable, logFile)

fprintf('Done.\n')
fprintf('Copied: %d\n', nCopied)
fprintf('Skipped/failed: %d\n', nSkipped)
fprintf('Log saved to:\n%s\n', logFile)


function cleanName = sanitize_folder_name(folderName)

cleanName = folderName;

badChars = {'<', '>', ':', '"', '/', '\', '|', '?', '*'};
for k = 1:numel(badChars)
    cleanName = strrep(cleanName, badChars{k}, '_');
end

cleanName = strtrim(cleanName);

if isempty(cleanName)
    cleanName = 'unnamed_recording';
end

end


function sessionName = get_session_folder_name(filePath, sourceRoot)

relativePath = strrep(filePath, sourceRoot, '');
if startsWith(relativePath, filesep)
    relativePath = relativePath(2:end);
end

parts = strsplit(relativePath, filesep);

if isempty(parts)
    sessionName = 'unknown_session';
else
    sessionName = parts{1};
end

sessionName = sanitize_folder_name(sessionName);

end


function finalFolder = make_unique_folder_if_needed(folderPath)

finalFolder = folderPath;

if ~exist(finalFolder, 'dir')
    return
end

k = 2;
while exist(finalFolder, 'dir')
    finalFolder = sprintf('%s__copy%d', folderPath, k);
    k = k + 1;
end

end