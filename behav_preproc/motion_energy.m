function [motion_tsd, t, motion_energy, roi_rect, info] = motion_energy(datapath)
% motion_energy
% Compute motion energy inside a fixed ROI for a video in datapath/video.
%
% INPUT
%   datapath : session folder (string), containing subfolder 'video'
%
% OUTPUT
%   motion_tsd    : tsd object with motion energy vs time
%   t             : [nFrames-1 x 1] time in seconds (starting at 0)
%   motion_energy : [nFrames-1 x 1] motion energy values
%   roi_rect      : [x y w h] ROI bounding box (pixels)
%   info          : struct with metadata (ROI, video, avg map, etc.)

% -------------------------------------------------------------------------
% 1) Find video file
% -------------------------------------------------------------------------
video_dir = fullfile(datapath, 'video');
file_list = dir(fullfile(video_dir, '*004*'));
if isempty(file_list)
    error('motion_energy:NoVideo', ...
        'No video file matching *004* found in %s', video_dir);
end

file = file_list(1); % take the first match
video_file = fullfile(file.folder, file.name);

v = VideoReader(video_file);

if ~hasFrame(v)
    error('motion_energy:EmptyVideo', 'Video has no frames: %s', video_file);
end

% -------------------------------------------------------------------------
% 2) Read first frame, make it double, get size
% -------------------------------------------------------------------------
frame1 = readFrame(v);
if ndims(frame1) == 3
    % assume already "black and white" encoded as 3 identical channels
    frame1 = frame1(:,:,1);
end
frame1 = double(frame1);

[H, W] = size(frame1);

% -------------------------------------------------------------------------
% 3) Load or create ROI mask (per animal, stored in parent folder)
% -------------------------------------------------------------------------
parent_dir = fileparts(datapath);
default_roi_file = fullfile(parent_dir, 'motion_energy_roi.mat');

if exist(default_roi_file, 'file')
    S = load(default_roi_file);
    if isfield(S, 'roi_struct')
        roi_struct = S.roi_struct;
    elseif isfield(S, 'roi')
        roi_struct = S.roi;
    else
        error('motion_energy:BadRoiFile', ...
            'ROI file %s does not contain roi_struct or roi.', default_roi_file);
    end
else
    % if no ROI file, make it (user draws ROI and chooses where to save)
    [roi_struct, default_roi_file] = make_motion_energy_roi(datapath);
end

if ~isfield(roi_struct, 'mask')
    error('motion_energy:NoMask', 'ROI struct has no field "mask".');
end

roi_mask_full = roi_struct.mask;
if ~isequal(size(roi_mask_full), [H W])
    error('motion_energy:MaskSizeMismatch', ...
        'ROI mask size [%d %d] does not match frame size [%d %d].', ...
        size(roi_mask_full,1), size(roi_mask_full,2), H, W);
end

% -------------------------------------------------------------------------
% 4) Get bounding box of ROI and cropped mask
% -------------------------------------------------------------------------
[rows, cols] = find(roi_mask_full);
if isempty(rows)
    error('motion_energy:EmptyRoi', 'ROI mask is empty.');
end

y1 = min(rows);
y2 = max(rows);
x1 = min(cols);
x2 = max(cols);

roi_rect = [x1, y1, x2 - x1 + 1, y2 - y1 + 1];
roi_mask = roi_mask_full(y1:y2, x1:x2);

% first ROI frame
prev_roi = frame1(y1:y2, x1:x2);

% -------------------------------------------------------------------------
% 5) Preallocate and loop through frames
% -------------------------------------------------------------------------
n_est = max(0, floor(v.Duration * v.FrameRate) - 1);
motion_energy = zeros(n_est, 1);
t = zeros(n_est, 1);
avg_map_sum = zeros(size(prev_roi));

idx = 0;

while hasFrame(v)
    frame = readFrame(v);
    if ndims(frame) == 3
        frame = frame(:,:,1);  % first channel only (no rgb2gray)
    end
    frame = double(frame);

    curr_roi = frame(y1:y2, x1:x2);

    diff_roi = abs(curr_roi - prev_roi);
    diff_roi = diff_roi .* roi_mask;  % zero outside mask

    idx = idx + 1;
    motion_energy(idx) = sum(diff_roi(:));
    t(idx) = v.CurrentTime;

    avg_map_sum = avg_map_sum + diff_roi;

    prev_roi = curr_roi;
end

% trim arrays
motion_energy = motion_energy(1:idx);
t = t(1:idx);

% shift time so it starts at 0
if ~isempty(t)
    t = t - t(1);
end

if idx > 0
    avg_map = avg_map_sum ./ idx;
else
    avg_map = avg_map_sum;
end

% -------------------------------------------------------------------------
% 6) Convert motion energy to tsd
% -------------------------------------------------------------------------
motion_tsd = tsd(t(:)*1e4, motion_energy(:));

% -------------------------------------------------------------------------
% 7) Build plots:
%   - time course
%   - power distribution (histogram)
%   - average motion energy map on ROI
% -------------------------------------------------------------------------
% Time course
figure;
subplot(3,1,1);
plot(Range(motion_tsd, 's'), Data(motion_tsd), '-');
xlabel('Time (s)');
ylabel('Motion energy');
title('Motion energy time course');

% Power distribution (histogram)
subplot(3,1,2);
histogram(motion_energy, 150);
xlabel('Motion energy');
ylabel('Count');
title('Motion energy distribution');

% Smoothed trace (optional, e.g. ~1 s window)
win = max(1, round(v.FrameRate));
smooth_me = movmean(motion_energy, win);

subplot(3,1,3);
plot(Range(motion_tsd, 's'), smooth_me, '-');
xlabel('Time (s)');
ylabel('Motion energy (smoothed)');
title(sprintf('Smoothed motion energy (window ~%g frames)', win));

% Average motion map on ROI
figure;
full_avg = zeros(H, W);
full_avg(y1:y2, x1:x2) = avg_map .* roi_mask;
imagesc(full_avg);
axis image off;
colorbar;
title('Average motion energy in ROI');

% -------------------------------------------------------------------------
% 8) Info struct
% -------------------------------------------------------------------------
info = struct();
info.video_file      = file.name;
info.video_folder    = file.folder;
info.frame_rate      = v.FrameRate;
info.roi_rect        = roi_rect;
info.roi_file        = default_roi_file;
info.roi_mask_full   = roi_mask_full;
info.avg_map_cropped = avg_map;
info.datapath        = datapath;
end
