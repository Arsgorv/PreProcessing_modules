function  [restored_data, shifts2, template2, options_nonrigid] = frames_alignment_AG(datapath, temp_data, plt, sess)
% Non-rigid motion correction of fUS frames using NoRMCorre.
% temp_data: 3D matrix [x y time] for a single slice
% plt:       0/1, show diagnostic plots and video
% sess:      index (used only in video name)

if nargin < 2 || isempty(plt)
    plt = 0;
end
if nargin < 3
    sess = 1;
end

% Step 1: NoRMCorre arguments
args = struct;
args.normcorreFilterType = 'defaultSVD';
args = normcorre_set_arguments_AG(args);

% Step 2: interactive cropping
bottom_adjust = 0;

while bottom_adjust <= 0
    figure;
    imagesc(temp_data(:,:,1));
    axis image; colormap hot;
    title('Draw ROI (imrect), double-click to confirm');

    roi = imrect;
    wait(roi);
    roiPosition = getPosition(roi);  % [x y width height]
    close;

    x_start  = round(roiPosition(1));
    y_start  = round(roiPosition(2));
    width    = round(roiPosition(3));
    height   = round(roiPosition(4));

    x_start = max(1, x_start);
    y_start = max(1, y_start);

    x_end = min(x_start + width  - 1, size(temp_data,2));
    y_end = min(y_start + height - 1, size(temp_data,1));

    crop_cols = x_start:x_end;
    crop_rows = y_start:y_end;

    left_adjust   = crop_cols(1) - 1;
    right_adjust  = size(temp_data,2) - crop_cols(end);
    top_adjust    = crop_rows(1) - 1;
    bottom_adjust = size(temp_data,1) - crop_rows(end);

    disp(['left margin:   ' num2str(left_adjust)]);
    disp(['right margin:  ' num2str(right_adjust)]);
    disp(['top margin:    ' num2str(top_adjust)]);
    disp(['bottom margin: ' num2str(bottom_adjust)]);

    if bottom_adjust <= 0
        disp('bottom_adjust <= 0, please redraw ROI');
    end
end

% Step 3: data for NoRMCorre
Y_in  = temp_data(crop_rows, crop_cols, :);
Y_out = Y_in;

args.normcorreArguments.grid_size = [size(Y_in,1), 6];

if plt
    figure;
    subplot(1,2,1); imagesc(Y_in(:,:,1));  axis image; title('Cropped');
    subplot(1,2,2); imagesc(temp_data(:,:,1)); axis image; title('Original');
end

% Step 4: run NoRMCorre
options_nonrigid = NoRMCorreSetParms( ...
    'd1',            size(Y_in,1), ...
    'd2',            size(Y_in,2), ...
    'grid_size',     args.normcorreArguments.grid_size, ...
    'mot_uf',        args.normcorreArguments.mot_uf, ...
    'max_shift',     args.normcorreArguments.max_shift, ...
    'max_dev',       args.normcorreArguments.max_dev, ...
    'init_batch',    args.normcorreArguments.init_batch, ...
    'overlap_pre',   args.normcorreArguments.overlap_pre, ...
    'overlap_post',  args.normcorreArguments.overlap_post, ...
    'shifts_method', args.normcorreArguments.shifts_method, ...
    'min_diff',      args.normcorreArguments.min_diff, ...
    'correct_bidir', args.normcorreArguments.correct_bidir);

[aligned_data, shifts2, template2, options_nonrigid] = ...
    normcorre_batch_twoArrays(Y_in, Y_out, options_nonrigid);

% Step 5: restore full frame
restored_data = temp_data;
restored_data(crop_rows, crop_cols, :) = aligned_data;

if plt
    figure;
    subplot(1,3,1); imagesc(temp_data(:,:,1));   axis image; title('Original');
    subplot(1,3,2); imagesc(aligned_data(:,:,1)); axis image; title('Aligned (crop)');
    subplot(1,3,3); imagesc(restored_data(:,:,1)); axis image; title('Restored');
end

% Step 6: optional diagnostic video
if plt
    nnY = quantile(Y_out(:), 0.1);
    mmY = quantile(Y_out(:), 0.9);

    [cM2,~,~] = motion_metrics(restored_data, 10);
    T = numel(cM2);

    demoFig = figure('visible', 'on');
    set(demoFig,'Position',[50 50 1200 800]);

    writerObj = VideoWriter([datapath filesep 'fUS' filesep 'video_alignment_session_' num2str(sess) '_PostExp']);
    set(writerObj,'FrameRate',10);
    open(writerObj);
    t_start = max(1, floor(size(restored_data,3) * 0.99));
    for t = t_start:size(restored_data,3)
        subplot(1,3,1);
        imagesc(sqrt(temp_data(:,:,t)), sqrt([nnY, mmY]));
        axis equal; axis tight;
        title(sprintf('Raw %i/%i', t, T));

        subplot(1,3,2);
        imagesc(sqrt(restored_data(:,:,t)), sqrt([nnY, mmY]));
        axis equal; axis tight;
        title('Corrected');
        colormap('hot');

        subplot(1,3,3);
        if t <= numel(shifts2)
            quiv = quiver( ...
                flipud(shifts2(t).shifts_up(:,:,1,2)), ...
                flipud(shifts2(t).shifts_up(:,:,1,1)), ...
                'AutoScale', 'off');
            ylim(quiv.Parent,[-10 10]);
            axis equal;
            title('Shifts');
        end
        drawnow;

        frameOut = getframe(demoFig);
        writeVideo(writerObj, frameOut);
    end

    close(writerObj);
end
end