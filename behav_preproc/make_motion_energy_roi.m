function [roi_struct, roi_file_out] = make_motion_energy_roi(datapath)
% make_motion_energy_roi
% Helper to create and save an ROI mask for motion_energy().
%
% - derives animal name from datapath
% - opens a representative frame
% - asks user to draw ROI
% - asks user where to save the ROI mask (.mat)
%
% OUTPUT
%   roi_struct   : struct with mask, rect, metadata
%   roi_file_out : full path to saved .mat file

% -------------------------------------------------------------------------
% 1) Find a video and read one frame
% -------------------------------------------------------------------------
video_dir = fullfile(datapath, 'video');
file_list = dir(fullfile(video_dir, '*004*'));
if isempty(file_list)
    error('make_motion_energy_roi:NoVideo', ...
        'No video file matching *004* found in %s', video_dir);
end

file = file_list(1);
video_file = fullfile(file.folder, file.name);
v = VideoReader(video_file);

if ~hasFrame(v)
    error('make_motion_energy_roi:EmptyVideo', ...
        'Video has no frames: %s', video_file);
end

frame1 = readFrame(v);
if ndims(frame1) == 3
    frame1 = frame1(:,:,1);
end
frame1 = double(frame1);

[H, W] = size(frame1); %#ok<NASGU>

% -------------------------------------------------------------------------
% 2) Let user draw ROI
% -------------------------------------------------------------------------
figure;
imagesc(frame1);
colormap gray;
axis image off;
title('Draw ROI (use mouse), double-click inside to confirm');

h = imrect;
rect = wait(h);            % [x y w h]
mask = createMask(h);      % logical mask same size as frame
close;

% -------------------------------------------------------------------------
% 3) Derive animal name from datapath
% -------------------------------------------------------------------------
[parent_dir, session_name] = fileparts(datapath);
[~, animal_name] = fileparts(parent_dir);

% -------------------------------------------------------------------------
% 4) Ask user where to save ROI mask
% -------------------------------------------------------------------------
default_name = 'motion_energy_roi.mat';
[roi_file, roi_folder] = uiputfile(default_name, ...
    'Save ROI mask (.mat)', fullfile(parent_dir, default_name));

if isequal(roi_file, 0)
    error('make_motion_energy_roi:UserCancelled', ...
        'User cancelled ROI saving.');
end

roi_file_out = fullfile(roi_folder, roi_file);

% -------------------------------------------------------------------------
% 5) Build roi_struct and save
% -------------------------------------------------------------------------
roi_struct = struct();
roi_struct.mask            = mask;
roi_struct.rect            = rect;
roi_struct.animal          = animal_name;
roi_struct.session         = session_name;
roi_struct.datapath_example= datapath;
roi_struct.video_file      = file.name;
roi_struct.created         = datestr(now);

save(roi_file_out, 'roi_struct');

fprintf('ROI saved to %s\n', roi_file_out);
end
