function draw_masks(datapath)
% draw_masks_session
% Inspect and draw anatomical masks for one fUS session.
%
% Expects per-slice files:
%   datapath\fUS\RP_data_<tail>_slice_A.mat (and/or B,C,D)
% where each file contains at least:
%   - cat_tsd (preferred) or tsd_raw
%   - optionally masks struct
%   - optionally legacy masks: Hpc, tissue, far_out, vessel
%
% Produces/updates:
%   - struct masks inside each RP_data_* file
%   - PNG overview: datapath\figures\anatomy\masks_<tail>_slice_X.png

if nargin < 1
    error('draw_masks_session:NoPath', 'Provide datapath for a single session.');
end

% tail = last part of the path, e.g. 20220517_2_n_S
tail = regexp(datapath, '[^\\\/]+$', 'match', 'once');

% animal = folder just before tail, e.g. Edel/Kosichka
parts = regexp(datapath, '[\\\/]', 'split');
if numel(parts) >= 2
    animal_str = parts{end-1};
else
    animal_str = 'animal?';
end
session_str = tail;

figdir = fullfile(datapath, 'figures', 'anatomy');
if ~exist(figdir, 'dir')
    mkdir(figdir);
end

% ROI definitions
roiLabels = { ...
    'Primary ACx', ...
    'Secondary ACx', ...
    'Hippocampus', ...
    'Thalamus', ...
    'Tissue', ...
    'Vessel', ...
    'Far_out'};

% Field names in masks struct
% Keep existing names for these: Hpc, tissue, far_out, vessel
roiFields = { ...
    'ACx_p', ...      % Primary ACx
    'ACx_s', ...      % Secondary ACx
    'Hpc', ...        % Hippocampus (legacy: Hpc)
    'Thal', ...       % Thalamus
    'tissue', ...     % Tissue (legacy: tissue)
    'vessel', ...     % Vessel (legacy: vessel)
    'far_out'};       % Far_out (legacy: far_out)

roiColors = lines(numel(roiLabels));
sliceLetters = 'ABCD';

for s = 1:length(sliceLetters)
    slice_label = sliceLetters(s);

    rpFile = fullfile(datapath, 'fUS', ['RP_data_' tail '_slice_' slice_label '.mat']);
    if ~exist(rpFile, 'file')
        continue;
    end

    fprintf('\n===== Session %s | Slice %s =====\n', tail, slice_label);

    S = load(rpFile);  % may contain cat_tsd, tsd_raw, masks, legacy masks

    % ------ choose tsd (aligned if available) ------
    if isfield(S, 'cat_tsd')
        tsd_use = S.cat_tsd;
    elseif isfield(S, 'tsd_raw')
        tsd_use = S.tsd_raw;
    else
        error('File %s has neither cat_tsd nor tsd_raw.', rpFile);
    end

    D  = Data(tsd_use.data);     % [T x (Nx*Ny)]
    Nx = tsd_use.Nx;
    Ny = tsd_use.Ny;
    T  = size(D,1);

    movie = reshape(D', Nx, Ny, T);  % [x y time]

    % ------ construct background image ------
    % mean across time, then log10 with eps for safety
%     bg = log10(mean(movie, 3) + eps);
    bg = mean(movie, 3) + eps;

    % percentile-based contrast
    bg_vec = bg(~isnan(bg));
    if ~isempty(bg_vec)
        low  = prctile(bg_vec, 1);
        high = prctile(bg_vec, 98);
        clim = [low high];
    else
        clim = [];
    end

    % ------ build masks struct, migrate legacy vars ------
    if isfield(S, 'masks')
        masks = S.masks;
    else
        masks = struct();
    end

    % migrate legacy mask variables if present and not already in struct
    legacyNames = {'Hpc','tissue','far_out','vessel'};
    for ln = 1:numel(legacyNames)
        fname = legacyNames{ln};
        if isfield(S, fname)
            if ~isfield(masks, fname) || isempty(masks.(fname))
                masks.(fname) = S.(fname);
            end
        end
    end

    % ------ plot existing masks overview ------
    fig_exist = figure('Color','w','Units','normalized','Position',[0.05 0.05 0.7 0.8]);
    if isempty(clim)
        imagesc(bg);
    else
        imagesc(bg, clim);
    end
    axis image off;
    colormap(viridis);
    hold on;

    for r = 1:numel(roiFields)
        fName = roiFields{r};
        if isfield(masks, fName) && ~isempty(masks.(fName))
            thisMask = masks.(fName);
            if isequal(size(thisMask), size(bg))
                contour(thisMask, [0.5 0.5], 'Color', roiColors(r,:), 'LineWidth', 1.5);
            end
        end
    end

    title(sprintf('%s | %s | slice %s (existing masks)', ...
        animal_str, session_str, slice_label), 'Interpreter','none');

    out_png = fullfile(figdir, sprintf('masks_%s_slice_%s.png', tail, slice_label));
    saveas(fig_exist, out_png);
    close(fig_exist);

    % ------ interactive editing / drawing ------
    anyChanged = false;

    for r = 1:numel(roiLabels)
        roiLabel = roiLabels{r};
        roiField = roiFields{r};

        % does this ROI already have a mask?
        hasExisting = isfield(masks, roiField) && ~isempty(masks.(roiField)) ...
                      && isequal(size(masks.(roiField)), size(bg));

        % if legacy name and ROIField differ (e.g. Hpc)
        % hasExisting already counts them via migration above

        if hasExisting
            % ask what to do with existing mask
            choiceExist = questdlg( ...
                sprintf('ROI %s (slice %s) already has a mask. What do you want to do?', roiLabel, slice_label), ...
                'Existing mask', ...
                'Keep','Redraw','Clear','Keep');

            if strcmp(choiceExist, 'Keep')
                continue;
            elseif strcmp(choiceExist, 'Clear')
                masks.(roiField) = [];
                anyChanged = true;
                continue;
            elseif strcmp(choiceExist, 'Redraw')
                % fall through to drawing
            else
                % dialog closed
                continue;
            end
        else
            % no existing mask: ask if draw
            choice = questdlg( ...
                sprintf('Draw new mask for %s (slice %s)?', roiLabel, slice_label), ...
                'Masks', ...
                'Yes','Skip','Skip');
            if strcmp(choice, 'Skip')
                continue;
            end
        end

        % ------ choose which ROIs to show as background ------
        haveAnyMask = false(1, numel(roiFields));
        for rr = 1:numel(roiFields)
            fName2 = roiFields{rr};
            if isfield(masks, fName2) && ~isempty(masks.(fName2)) ...
                    && isequal(size(masks.(fName2)), size(bg))
                haveAnyMask(rr) = true;
            end
        end

        visibleIdx = find(haveAnyMask);
        if ~isempty(visibleIdx)
            labelList = roiLabels;
            defIdx    = visibleIdx;
            [selIdx, ok] = listdlg( ...
                'PromptString', sprintf('Select ROIs to show as background (slice %s):', slice_label), ...
                'ListString', labelList, ...
                'InitialValue', defIdx, ...
                'SelectionMode','multiple', ...
                'ListSize', [220 150]);
            if ok == 1 && ~isempty(selIdx)
                visibleIdx = selIdx(:)';
            end
        end

        % ------ draw until satisfied ------
        satisfied = false;
        while ~satisfied
            fig_draw = figure('Color','w','Units','normalized','Position',[0.05 0.05 0.7 0.8]);
            if isempty(clim)
                imagesc(bg);
            else
                imagesc(bg, clim);
            end
            axis image off; colormap(gray);
            hold on;

            % show selected background ROIs
            for rr = visibleIdx
                fName2 = roiFields{rr};
                thisMask2 = masks.(fName2);
                if isequal(size(thisMask2), size(bg))
                    contour(thisMask2, [0.5 0.5], 'Color', roiColors(rr,:), 'LineWidth', 1);
                end
            end

            title(sprintf('%s | %s | slice %s | ROI: %s', ...
                animal_str, session_str, slice_label, roiLabel), ...
                'Interpreter','none');

            disp(['Draw polygon for ' roiLabel ' (slice ' slice_label '), double-click to close.']);
            newMask = roipoly;

            % preview
            clf(fig_draw);
            if isempty(clim)
                imagesc(bg);
            else
                imagesc(bg, clim);
            end
            axis image off; colormap(gray);
            hold on;
            contour(newMask, [0.5 0.5], 'r', 'LineWidth', 1.5);
            title(sprintf('%s | %s | slice %s | ROI: %s (preview)', ...
                animal_str, session_str, slice_label, roiLabel), ...
                'Interpreter','none');

            choice2 = questdlg('Satisfied with this mask?', ...
                ['ROI: ' roiLabel], ...
                'Yes','Redraw','Redraw');

            close(fig_draw);

            if strcmp(choice2, 'Yes')
                masks.(roiField) = newMask;
                anyChanged = true;
                satisfied = true;
            else
                satisfied = false;
            end
        end
    end

    % ------ save + replot combined if changed ------
    if anyChanged
        save(rpFile, 'masks', '-append');

        fig_all = figure('Color','w','Units','normalized','Position',[0.05 0.05 0.7 0.8]);
        if isempty(clim)
            imagesc(bg);
        else
            imagesc(bg, clim);
        end
        axis image off; colormap(gray);
        hold on;

        for r = 1:numel(roiFields)
            fName = roiFields{r};
            if isfield(masks, fName) && ~isempty(masks.(fName))
                thisMask = masks.(fName);
                if isequal(size(thisMask), size(bg))
                    contour(thisMask, [0.5 0.5], 'Color', roiColors(r,:), 'LineWidth', 1.5);
                end
            end
        end

        title(sprintf('%s | %s | slice %s (updated masks)', ...
            animal_str, session_str, slice_label), 'Interpreter','none');

        out_png = fullfile(figdir, sprintf('masks_%s_slice_%s.png', tail, slice_label));
        saveas(fig_all, out_png);
        close(fig_all);
    end
end
end