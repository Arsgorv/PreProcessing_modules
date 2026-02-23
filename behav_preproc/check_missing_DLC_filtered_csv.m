function missing = check_missing_DLC_filtered_csv(rootDir)
% check_missing_DLC_filtered_csv
% Checks presence of DLC *_filtered.csv for FACE (24934004) and EYE (24934007)
% in each session's "video" folder under rootDir.
%
% OUTPUT
%   missing : table with columns {session_path, missing_face, missing_eye}
%
% Example:
%   rootDir = 'Z:\Arsenii\React_Active\training\Tvorozhok';
%   missing = RA_check_missing_filtered_csv(rootDir);

if nargin < 1 || isempty(rootDir) || ~ischar(rootDir)
    error('Provide rootDir as a char path.');
end
if ~exist(rootDir,'dir')
    error('Folder not found: %s', rootDir);
end

sess = list_session_dirs(rootDir);

session_path = {};
missing_face = [];
missing_eye  = [];

for i = 1:numel(sess)
    datapath = sess{i};
    videoDir = fullfile(datapath,'video');
    if ~exist(videoDir,'dir')
        continue
    end

    % list all filtered csvs in video folder
    ff = dir(fullfile(videoDir,'*_filtered.csv'));
    names = {ff.name};

    hasFace = any(contains(names,'24934004') & contains(names,'DLC'));
    hasEye  = any(contains(names,'24934007') & contains(names,'DLC'));

    if ~hasFace || ~hasEye
        session_path{end+1,1} = datapath; %#ok<AGROW>
        missing_face(end+1,1) = ~hasFace; %#ok<AGROW>
        missing_eye(end+1,1)  = ~hasEye;  %#ok<AGROW>
    end
end

missing = table(session_path, missing_face, missing_eye);

% Print concise report
if isempty(missing)
    fprintf('All sessions under %s have FACE and EYE *_filtered.csv\n', rootDir);
else
    fprintf('Missing *_filtered.csv in %d session(s):\n', height(missing));
    for i = 1:height(missing)
        cams = {};
        if missing.missing_face(i), cams{end+1} = 'FACE(24934004)'; end %#ok<AGROW>
        if missing.missing_eye(i),  cams{end+1} = 'EYE(24934007)';  end %#ok<AGROW>
        fprintf('  %s  -> missing: %s\n', missing.session_path{i}, strjoin(cams, ', '));
    end
end

end

% ---------------- helpers ----------------

function sess = list_session_dirs(rootDir)
% Recursively find session folders that contain a "video" subfolder.
sess = {};

d = dir(rootDir);
d = d([d.isdir]);
names = {d.name};
keep = ~ismember(names,{'.','..'});
d = d(keep);

for i = 1:numel(d)
    p = fullfile(rootDir, d(i).name);

    if exist(fullfile(p,'video'),'dir')
        sess{end+1,1} = p; %#ok<AGROW>
        continue
    end

    % recurse
    sub = list_session_dirs(p);
    if ~isempty(sub)
        sess = [sess; sub]; %#ok<AGROW>
    end
end
end
