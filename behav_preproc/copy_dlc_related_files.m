function copy_dlc_related_files(targetRoot, sourceRoot, overwrite)
% copy_dlc_related_files.m
%
% Copy DLC-related files from sourceRoot into each video's folder under targetRoot.
% Matching rule: any file in sourceRoot whose filename starts with the video base name.
% Example:
%   video base name: Mochi_Basler_acA800-510um__24934007__20251125_201351184
%   matches:         Mochi_Basler_acA800-510um__24934007__20251125_201351184DLC_*.csv/.h5/.pickle/... etc
%
% MATLAB R2018b compatible (uses dir with ** for recursion).

if nargin < 3 || isempty(overwrite)
    overwrite = false;
end

assert(ischar(targetRoot) || isstring(targetRoot), 'targetRoot must be a path string');
assert(ischar(sourceRoot) || isstring(sourceRoot), 'sourceRoot must be a path string');

targetRoot = char(targetRoot);
sourceRoot = char(sourceRoot);

if ~exist(targetRoot, 'dir')
    error('targetRoot does not exist: %s', targetRoot);
end
if ~exist(sourceRoot, 'dir')
    error('sourceRoot does not exist: %s', sourceRoot);
end

% ---- collect all videos under targetRoot (recursive) ----
videoExt = {'.mp4','.avi','.mov','.mkv','.mpg','.mpeg','.wmv','.m4v'};

vidFiles = [];
for i = 1:numel(videoExt)
    d = dir(fullfile(targetRoot, '**', ['*' videoExt{i}]));
    if ~isempty(d)
        vidFiles = [vidFiles; d]; %#ok<AGROW>
    end
end

if isempty(vidFiles)
    fprintf('No video files found under: %s\n', targetRoot);
    return;
end

fprintf('Found %d video(s) under targetRoot.\n', numel(vidFiles));

nCopiedTotal = 0;
nMissingTotal = 0;

for v = 1:numel(vidFiles)

    videoPath = fullfile(vidFiles(v).folder, vidFiles(v).name);
    [videoDir, baseName, ext] = fileparts(videoPath); %#ok<ASGLU>

    % Find all files in sourceRoot starting with baseName
    matches = dir(fullfile(sourceRoot, '**', [baseName '*']));

    % Remove directories
    matches = matches(~[matches.isdir]);

    % Optionally skip copying the raw video itself (if present in sourceRoot)
    % i.e. exact same name (baseName + a video extension)
    skipIdx = false(size(matches));
    for k = 1:numel(matches)
        [~, mBase, mExt] = fileparts(matches(k).name);
        if strcmp(mBase, baseName) && any(strcmpi(mExt, videoExt))
            skipIdx(k) = true;
        end
    end
    matches = matches(~skipIdx);

    if isempty(matches)
        fprintf('[MISS] %s  -> no related files found in sourceRoot\n', baseName);
        nMissingTotal = nMissingTotal + 1;
        continue;
    end

    fprintf('[%d/%d] %s  -> %d related file(s)\n', v, numel(vidFiles), baseName, numel(matches));

    for k = 1:numel(matches)
        srcPath = fullfile(matches(k).folder, matches(k).name);
        dstPath = fullfile(videoDir, matches(k).name);

        if exist(dstPath, 'file') && ~overwrite
            % If same size, skip quietly; if different size, still skip (safer) unless overwrite=true
            dstInfo = dir(dstPath);
            if ~isempty(dstInfo) && dstInfo(1).bytes == matches(k).bytes
                continue;
            else
                fprintf('  [SKIP exists] %s\n', matches(k).name);
                continue;
            end
        end

        % Ensure destination folder exists (should, but keep robust)
        if ~exist(videoDir, 'dir')
            mkdir(videoDir);
        end

        copyfile(srcPath, dstPath);
        nCopiedTotal = nCopiedTotal + 1;
    end
end

fprintf('\nDone.\n');
fprintf('Videos with no matches: %d\n', nMissingTotal);
fprintf('Total files copied:     %d\n', nCopiedTotal);

end
