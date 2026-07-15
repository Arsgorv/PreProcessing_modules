function plan = restructure_ephys_old_to_new_layout(datapath, varargin)

p = inputParser;

addParameter(p, 'execute', false);
addParameter(p, 'sourceFolders', {'FMA_OB', 'Neuropixels'});
addParameter(p, 'phaseNames', {'PreSleep', 'Conditioning', 'PostSleep', 'PostTest'});
addParameter(p, 'taskTag', 'RA');
addParameter(p, 'animalName', '');
addParameter(p, 'dateString', '');
addParameter(p, 'targetRoot', '');
addParameter(p, 'maxSearchDepth', 8);
addParameter(p, 'allowExistingTargets', false);
addParameter(p, 'checkSizesAfterCopy', true);
addParameter(p, 'verbose', true);

parse(p, varargin{:});
opts = p.Results;

opts.sourceFolders = make_cellstr(opts.sourceFolders);
opts.phaseNames = make_cellstr(opts.phaseNames);

ephysRoots = get_ephys_roots_from_input(datapath);

plan = struct([]);

for iRoot = 1:length(ephysRoots)
    thisRoot = ephysRoots{iRoot};
    
    thisPlan = build_restructure_plan(thisRoot, opts);
    thisPlan = validate_restructure_plan(thisPlan, opts);
    
    if opts.verbose
        print_plan_summary(thisPlan);
    end
    
    if opts.execute
        if ~thisPlan.valid
            error('Plan is not valid. Nothing was copied. Inspect plan.issues.');
        end
        
        thisPlan = execute_restructure_plan(thisPlan, opts);
    end
    
    plan = [plan thisPlan];
end

end


function plan = build_restructure_plan(ephysRoot, opts)

ephysRoot = remove_trailing_filesep(ephysRoot);

if isempty(opts.targetRoot)
    targetRoot = ephysRoot;
else
    targetRoot = remove_trailing_filesep(opts.targetRoot);
end

[defaultAnimal, defaultDate] = get_animal_and_date_from_path(ephysRoot);

if isempty(opts.animalName)
    animalName = defaultAnimal;
else
    animalName = opts.animalName;
end

if isempty(opts.dateString)
    dateString = defaultDate;
else
    dateString = normalize_date_string(opts.dateString);
end

sourceFolders = get_existing_source_folders(ephysRoot, opts.sourceFolders);
phaseNames = detect_phase_names(ephysRoot, sourceFolders, opts.phaseNames);

plan = struct();
plan.ephysRoot = ephysRoot;
plan.targetRoot = targetRoot;
plan.execute = opts.execute;
plan.valid = false;
plan.issues = {};
plan.operations = make_empty_operation_struct();
plan.sourceFolders = sourceFolders;
plan.phaseNames = phaseNames;
plan.animalName = animalName;
plan.dateString = dateString;
plan.taskTag = opts.taskTag;

opCount = 0;

for iPhase = 1:length(phaseNames)
    phaseName = phaseNames{iPhase};
    
    phaseRecords = struct([]);
    phaseSessionNames = {};
    
    for iSource = 1:length(sourceFolders)
        sourceFolder = sourceFolders{iSource};
        phaseSourceRoot = fullfile(ephysRoot, sourceFolder, phaseName);
        
        if ~exist(phaseSourceRoot, 'dir')
            continue
        end
        
        sessionDirs = list_dirs(phaseSourceRoot);
        
        if isempty(sessionDirs)
            contFolders = find_continuous_folders(phaseSourceRoot, 0, opts.maxSearchDepth);
            
            for iCont = 1:length(contFolders)
                rec = struct();
                rec.sourceFolder = sourceFolder;
                rec.sourceSession = phaseName;
                rec.continuousPath = contFolders{iCont};
                phaseRecords = [phaseRecords rec];
                phaseSessionNames{end + 1} = phaseName;
            end
        else
            for iSess = 1:length(sessionDirs)
                sessionPath = fullfile(phaseSourceRoot, sessionDirs(iSess).name);
                contFolders = find_continuous_folders(sessionPath, 0, opts.maxSearchDepth);
                
                for iCont = 1:length(contFolders)
                    rec = struct();
                    rec.sourceFolder = sourceFolder;
                    rec.sourceSession = sessionDirs(iSess).name;
                    rec.continuousPath = contFolders{iCont};
                    phaseRecords = [phaseRecords rec];
                    phaseSessionNames{end + 1} = sessionDirs(iSess).name;
                end
            end
        end
    end
    
    if isempty(phaseRecords)
        continue
    end
    
    targetSessionName = make_target_session_name(animalName, dateString, opts.taskTag, phaseName, phaseSessionNames);
    targetContinuousPath = fullfile(targetRoot, targetSessionName, 'recording1', 'continuous');
    
    for iRec = 1:length(phaseRecords)
        items = list_items(phaseRecords(iRec).continuousPath);
        
        for iItem = 1:length(items)
            src = fullfile(phaseRecords(iRec).continuousPath, items(iItem).name);
            dst = fullfile(targetContinuousPath, items(iItem).name);
            
            opCount = opCount + 1;
            op = struct();
            op.phase = phaseName;
            op.sourceFolder = phaseRecords(iRec).sourceFolder;
            op.sourceSession = phaseRecords(iRec).sourceSession;
            op.sourceContinuousPath = phaseRecords(iRec).continuousPath;
            op.itemName = items(iItem).name;
            op.sourcePath = src;
            op.destinationPath = dst;
            op.targetSessionName = targetSessionName;
            op.status = 'planned';
            op.message = '';
            
            plan.operations(opCount) = op;
        end
    end
end

end


function plan = validate_restructure_plan(plan, opts)

issues = {};

if ~exist(plan.ephysRoot, 'dir')
    issues{end + 1} = ['Ephys root does not exist: ' plan.ephysRoot];
end

if isempty(plan.sourceFolders)
    issues{end + 1} = 'No old source folders found.';
end

if isempty(plan.operations)
    issues{end + 1} = 'No copy operations were found. Check source folder names, phase names, and maxSearchDepth.';
end

destinationList = {};

for iOp = 1:length(plan.operations)
    src = plan.operations(iOp).sourcePath;
    dst = plan.operations(iOp).destinationPath;
    
    if ~exist(src, 'file') && ~exist(src, 'dir')
        issues{end + 1} = ['Missing source: ' src];
    end
    
    if path_is_inside(dst, src)
        issues{end + 1} = ['Unsafe destination inside source: ' dst];
    end
    
    if exist(dst, 'file') || exist(dst, 'dir')
        if opts.allowExistingTargets
            plan.operations(iOp).status = 'skip_existing';
            plan.operations(iOp).message = 'Destination already exists.';
        else
            issues{end + 1} = ['Destination already exists: ' dst];
        end
    end
    
    destinationList{end + 1} = lower(dst);
end

if ~isempty(destinationList)
    [uniqueDest, ~, idx] = unique(destinationList);
    
    for iDest = 1:length(uniqueDest)
        n = sum(idx == iDest);
        
        if n > 1
            repeated = find(idx == iDest);
            firstOp = repeated(1);
            issues{end + 1} = ['Duplicate destination in plan: ' plan.operations(firstOp).destinationPath];
        end
    end
end

plan.issues = issues;
plan.valid = isempty(issues);

end


function plan = execute_restructure_plan(plan, opts)

logFolder = fullfile(plan.targetRoot, 'logs');

if ~exist(logFolder, 'dir')
    mkdir(logFolder);
end

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
logFile = fullfile(logFolder, ['restructure_ephys_' timestamp '.txt']);

write_plan_log(plan, logFile, 'before_copy');

for iOp = 1:length(plan.operations)
    src = plan.operations(iOp).sourcePath;
    dst = plan.operations(iOp).destinationPath;
    
    if strcmp(plan.operations(iOp).status, 'skip_existing')
        continue
    end
    
    if exist(dst, 'file') || exist(dst, 'dir')
        error(['Destination appeared during execution. Aborting before overwrite: ' dst]);
    end
    
    dstParent = fileparts(dst);
    
    if ~exist(dstParent, 'dir')
        mkdir(dstParent);
    end
    
    fprintf('Copying %d/%d:\n%s\n--> %s\n', iOp, length(plan.operations), src, dst);
    
    [ok, msg] = copyfile(src, dst);
    
    if ~ok
        plan.operations(iOp).status = 'copy_failed';
        plan.operations(iOp).message = msg;
        write_plan_log(plan, logFile, 'copy_failed');
        error(['Copy failed: ' msg]);
    end
    
    if opts.checkSizesAfterCopy
        srcBytes = get_path_size_bytes(src);
        dstBytes = get_path_size_bytes(dst);
        
        if srcBytes ~= dstBytes
            plan.operations(iOp).status = 'size_mismatch';
            plan.operations(iOp).message = ['Source bytes: ' num2str(srcBytes) ', destination bytes: ' num2str(dstBytes)];
            write_plan_log(plan, logFile, 'size_mismatch');
            error(['Size mismatch after copy: ' dst]);
        end
    end
    
    plan.operations(iOp).status = 'copied';
    plan.operations(iOp).message = 'OK';
end

write_plan_log(plan, logFile, 'after_copy');

matLogFile = fullfile(logFolder, ['restructure_ephys_' timestamp '.mat']);
save(matLogFile, 'plan');

fprintf('\nFinished safely.\n');
fprintf('Log file: %s\n', logFile);
fprintf('MAT plan: %s\n', matLogFile);

end


function ephysRoots = get_ephys_roots_from_input(datapath)

ephysRoots = {};

if ischar(datapath)
    rootPath = normalize_ephys_root(datapath);
    ephysRoots{end + 1} = rootPath;
    return
end

if isstruct(datapath)
    for i = 1:numel(datapath)
        candidatePaths = get_candidate_paths_from_struct(datapath(i));
        
        for j = 1:length(candidatePaths)
            rootPath = normalize_ephys_root(candidatePaths{j});
            
            if exist(rootPath, 'dir')
                ephysRoots{end + 1} = rootPath;
            end
        end
    end
end

ephysRoots = unique(ephysRoots);

if isempty(ephysRoots)
    error('Could not find an ephys folder from datapath input.');
end

end


function candidatePaths = get_candidate_paths_from_struct(s)

candidatePaths = {};
fieldsToTry = {'ephys', 'Ephys', 'ephysPath', 'EphysPath', 'ephys_path', ...
    'path_ephys', 'ephysDir', 'EphysDir', 'datapath', 'DataPath', ...
    'path', 'Path', 'folder', 'Folder', 'basepath', 'BasePath'};

for iField = 1:length(fieldsToTry)
    fld = fieldsToTry{iField};
    
    if isfield(s, fld)
        val = s.(fld);
        
        if ischar(val)
            candidatePaths{end + 1} = val;
        elseif iscell(val)
            for k = 1:numel(val)
                if ischar(val{k})
                    candidatePaths{end + 1} = val{k};
                end
            end
        end
    end
end

end


function rootPath = normalize_ephys_root(pathIn)

pathIn = remove_trailing_filesep(pathIn);

if exist(pathIn, 'file') && ~exist(pathIn, 'dir')
    pathIn = fileparts(pathIn);
end

[~, lastName] = fileparts(pathIn);

if strcmpi(lastName, 'ephys')
    rootPath = pathIn;
    return
end

if exist(fullfile(pathIn, 'ephys'), 'dir')
    rootPath = fullfile(pathIn, 'ephys');
    return
end

rootPath = pathIn;

end


function sourceFolders = get_existing_source_folders(ephysRoot, requestedSourceFolders)

sourceFolders = {};

for i = 1:length(requestedSourceFolders)
    thisFolder = fullfile(ephysRoot, requestedSourceFolders{i});
    
    if exist(thisFolder, 'dir')
        sourceFolders{end + 1} = requestedSourceFolders{i};
    end
end

if ~isempty(sourceFolders)
    return
end

rootDirs = list_dirs(ephysRoot);

for iDir = 1:length(rootDirs)
    dirName = rootDirs(iDir).name;
    
    if strcmpi(dirName, 'logs')
        continue
    end
    
    testPath = fullfile(ephysRoot, dirName);
    childDirs = list_dirs(testPath);
    
    if ~isempty(childDirs)
        sourceFolders{end + 1} = dirName;
    end
end

end


function phaseNames = detect_phase_names(ephysRoot, sourceFolders, requestedPhaseNames)

phaseNames = requestedPhaseNames;

for iSource = 1:length(sourceFolders)
    sourcePath = fullfile(ephysRoot, sourceFolders{iSource});
    
    if ~exist(sourcePath, 'dir')
        continue
    end
    
    d = list_dirs(sourcePath);
    
    for i = 1:length(d)
        if ~cell_contains(phaseNames, d(i).name)
            phaseNames{end + 1} = d(i).name;
        end
    end
end

end


function contFolders = find_continuous_folders(startPath, currentDepth, maxDepth)

contFolders = {};

if currentDepth > maxDepth
    return
end

[~, name] = fileparts(startPath);

if strcmpi(name, 'continuous')
    contFolders{end + 1} = startPath;
    return
end

d = list_dirs(startPath);

for i = 1:length(d)
    childPath = fullfile(startPath, d(i).name);
    childCont = find_continuous_folders(childPath, currentDepth + 1, maxDepth);
    
    for j = 1:length(childCont)
        contFolders{end + 1} = childCont{j};
    end
end

end


function targetName = make_target_session_name(animalName, dateString, taskTag, phaseName, phaseSessionNames)

[sessionDate, sessionTime] = get_earliest_session_datetime(phaseSessionNames);

if isempty(dateString)
    dateString = sessionDate;
end

parts = {};

if ~isempty(animalName)
    parts{end + 1} = animalName;
end

if ~isempty(dateString)
    parts{end + 1} = dateString;
end

if ~isempty(sessionTime)
    parts{end + 1} = sessionTime;
end

if ~isempty(taskTag)
    parts{end + 1} = taskTag;
end

parts{end + 1} = phaseName;

targetName = parts{1};

for i = 2:length(parts)
    targetName = [targetName '_' parts{i}];
end

targetName = sanitize_folder_name(targetName);

end


function [bestDate, bestTime] = get_earliest_session_datetime(sessionNames)

bestDate = '';
bestTime = '';
bestNum = inf;

for i = 1:length(sessionNames)
    name = sessionNames{i};
    tok = regexp(name, '(\d{4}-\d{2}-\d{2})[^0-9]+(\d{2})-(\d{2})-(\d{2})', 'tokens', 'once');
    
    if isempty(tok)
        continue
    end
    
    thisDate = tok{1};
    thisTimeForFolder = [tok{2} '-' tok{3} '-' tok{4}];
    thisTimeForDatenum = [tok{2} ':' tok{3} ':' tok{4}];
    
    try
        thisNum = datenum([thisDate ' ' thisTimeForDatenum], 'yyyy-mm-dd HH:MM:SS');
    catch
        thisNum = inf;
    end
    
    if thisNum < bestNum
        bestNum = thisNum;
        bestDate = thisDate;
        bestTime = thisTimeForFolder;
    end
end

end


function [animalName, dateString] = get_animal_and_date_from_path(ephysRoot)

animalName = '';
dateString = '';

parent1 = fileparts(remove_trailing_filesep(ephysRoot));
[parent2, dateFolder] = fileparts(parent1);
[~, animalFolder] = fileparts(parent2);

if ~isempty(animalFolder)
    animalName = animalFolder;
end

dateString = normalize_date_string(dateFolder);

end


function dateString = normalize_date_string(dateIn)

dateString = '';

if isempty(dateIn)
    return
end

if length(dateIn) == 8 && all(dateIn >= '0' & dateIn <= '9')
    dateString = [dateIn(1:4) '-' dateIn(5:6) '-' dateIn(7:8)];
    return
end

tok = regexp(dateIn, '(\d{4})-(\d{2})-(\d{2})', 'tokens', 'once');

if ~isempty(tok)
    dateString = [tok{1} '-' tok{2} '-' tok{3}];
end

end


function d = list_dirs(pathIn)

d0 = dir(pathIn);
d = d0([d0.isdir]);

keep = true(size(d));

for i = 1:length(d)
    if strcmp(d(i).name, '.') || strcmp(d(i).name, '..')
        keep(i) = false;
    end
end

d = d(keep);

end


function items = list_items(pathIn)

items0 = dir(pathIn);
keep = true(size(items0));

for i = 1:length(items0)
    if strcmp(items0(i).name, '.') || strcmp(items0(i).name, '..')
        keep(i) = false;
    end
end

items = items0(keep);

end


function bytes = get_path_size_bytes(pathIn)

bytes = 0;

if exist(pathIn, 'file') && ~exist(pathIn, 'dir')
    d = dir(pathIn);
    
    if ~isempty(d)
        bytes = d.bytes;
    end
    
    return
end

if exist(pathIn, 'dir')
    items = list_items(pathIn);
    
    for i = 1:length(items)
        childPath = fullfile(pathIn, items(i).name);
        
        if items(i).isdir
            bytes = bytes + get_path_size_bytes(childPath);
        else
            bytes = bytes + items(i).bytes;
        end
    end
end

end


function tf = path_is_inside(childPath, parentPath)

childPath = lower(remove_trailing_filesep(childPath));
parentPath = lower(remove_trailing_filesep(parentPath));

if length(childPath) <= length(parentPath)
    tf = false;
    return
end

tf = strcmp(childPath(1:length(parentPath)), parentPath);

end


function out = remove_trailing_filesep(in)

out = in;

while ~isempty(out)
    lastChar = out(end);
    
    if lastChar == filesep || lastChar == '/' || lastChar == '\'
        out(end) = [];
    else
        break
    end
end

end


function c = make_cellstr(x)

if ischar(x)
    c = {x};
else
    c = x;
end

end


function tf = cell_contains(c, value)

tf = false;

for i = 1:length(c)
    if strcmpi(c{i}, value)
        tf = true;
        return
    end
end

end


function nameOut = sanitize_folder_name(nameIn)

nameOut = nameIn;
badChars = '<>:"/\|?*';

for i = 1:length(badChars)
    nameOut(nameOut == badChars(i)) = '_';
end

end


function print_plan_summary(plan)

fprintf('\nEphys root:\n%s\n', plan.ephysRoot);
fprintf('\nTarget root:\n%s\n', plan.targetRoot);
fprintf('\nSource folders:\n');

for i = 1:length(plan.sourceFolders)
    fprintf('  %s\n', plan.sourceFolders{i});
end

fprintf('\nPhases:\n');

for i = 1:length(plan.phaseNames)
    fprintf('  %s\n', plan.phaseNames{i});
end

fprintf('\nCopy operations planned: %d\n', length(plan.operations));

if plan.valid
    fprintf('Plan status: VALID\n');
else
    fprintf('Plan status: INVALID\n');
    fprintf('\nIssues:\n');
    
    for i = 1:length(plan.issues)
        fprintf('  %s\n', plan.issues{i});
    end
end

fprintf('\nDry run only unless execute=true.\n');

end


function write_plan_log(plan, logFile, stage)

fid = fopen(logFile, 'a');

if fid < 0
    warning(['Could not write log file: ' logFile]);
    return
end

fprintf(fid, '\n===== %s =====\n', stage);
fprintf(fid, 'Time: %s\n', datestr(now));
fprintf(fid, 'Ephys root: %s\n', plan.ephysRoot);
fprintf(fid, 'Target root: %s\n', plan.targetRoot);
fprintf(fid, 'Valid: %d\n', plan.valid);

if ~isempty(plan.issues)
    fprintf(fid, '\nIssues:\n');
    
    for i = 1:length(plan.issues)
        fprintf(fid, '%s\n', plan.issues{i});
    end
end

fprintf(fid, '\nOperations:\n');

for i = 1:length(plan.operations)
    op = plan.operations(i);
    fprintf(fid, '\nOperation %d\n', i);
    fprintf(fid, 'Phase: %s\n', op.phase);
    fprintf(fid, 'Source folder: %s\n', op.sourceFolder);
    fprintf(fid, 'Source session: %s\n', op.sourceSession);
    fprintf(fid, 'Item: %s\n', op.itemName);
    fprintf(fid, 'Source: %s\n', op.sourcePath);
    fprintf(fid, 'Destination: %s\n', op.destinationPath);
    fprintf(fid, 'Status: %s\n', op.status);
    fprintf(fid, 'Message: %s\n', op.message);
end

fclose(fid);

end

function operations = make_empty_operation_struct()

operations = struct( ...
    'phase', {}, ...
    'sourceFolder', {}, ...
    'sourceSession', {}, ...
    'sourceContinuousPath', {}, ...
    'itemName', {}, ...
    'sourcePath', {}, ...
    'destinationPath', {}, ...
    'targetSessionName', {}, ...
    'status', {}, ...
    'message', {});

end