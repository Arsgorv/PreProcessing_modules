function cropping_prior_to_alignment(datapath)

tail = regexp(datapath, '[^\\\/]+$', 'match', 'once');

sliceLetters = 'ABCD';

for s = 1:length(sliceLetters)
    slice_label = sliceLetters(s);
    rpFile = [datapath filesep 'fUS' filesep 'RP_data_' tail '_slice_' slice_label '.mat'];

    if ~exist(rpFile, 'file')
        % No such slice in this session, skip
        continue;
    end

    disp(['Processing ' rpFile])

    % Expect tsd_raw inside
    S = load(rpFile);
    if ~isfield(S, 'tsd_raw')
        warning('File %s does not contain tsd_raw. Skipping.', rpFile);
        continue;
    end
    tsd_raw = S.tsd_raw;

    % tsd_raw: time x (Nx*Ny)
    D  = Data(tsd_raw.data);      % [T x (Nx*Ny)]
    Nx = tsd_raw.Nx;
    Ny = tsd_raw.Ny;
    T  = size(D,1);

    % movie: [x y time]
    temp_data = reshape(D', Nx, Ny, T);

    % interactive cropping
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

    % append to RP_<tail>_X.mat
    save(rpFile, 'crop_cols', 'crop_rows', '-append');
end
end