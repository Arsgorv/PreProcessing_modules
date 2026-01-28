function behaviour_preprocessing(datapath)
% This script computes basic parameters and generates figures from DLC tracking data.
% It uses synchronized DLC data (with the time column aligned to openephys triggers)
% so that the ephys (Brain/LFP) and video/DLC signals are on the same time base.

%% Initialization
[~, Session_params.session_selection, ~] = fileparts(datapath);

Session_params.do_plots = 0;          % set 0 to skip QC plots
Session_params.save_plots = 0;        % save pngs in video/
Session_params.qc_subsample = 20;     % for scatter decimation
Session_params.fig_visibility = 'on';

Session_params.animal_name = 'Unknown';
animals = {'Shropshire','Brynza','Labneh','Tvorozhok','Kosichka','Mochi','Edel','Chabichou','Ficello','Kiri','Brayon'};
for i=1:numel(animals)
    if contains(datapath, animals{i})
        Session_params.animal_name = animals{i};
        break
    end
end

disp(['Processing session: '  Session_params.animal_name ' ' Session_params.session_selection]);

%% Load DLC Synchronized Data
dlc_file{1} = dir(fullfile(datapath, 'video', 'synchronized_DLC_data_face.csv'));
dlc_file{2} = dir(fullfile(datapath, 'video', 'synchronized_DLC_data_eye.csv'));
if isempty(dlc_file{1}) || isempty(dlc_file{2})
    error('One or both DLC files not found in %s/video', datapath);
end

filename{1} = dlc_file{1}.name;
filename{2} = dlc_file{2}.name;
disp(['DLC data face: ' filename{1} ' ; DLC data eye:' filename{2}])

dataStruct{1} = importdata(fullfile(datapath, 'video', filename{1}));
dataStruct{2} = importdata(fullfile(datapath, 'video', filename{2}));
for k = 1:2
    if isstruct(dataStruct{k})
        data{k} = dataStruct{k}.data;
    else
        data{k} = dataStruct{k};
    end
end

% Verify that the first column (time) is monotonic
for id = 1:2
    time_vals = data{id}(:,1);
    if any(diff(time_vals) < 0)
        warning('Time stamps in file %d are not monotonic. Sorting...', id);
        [~, idx] = sort(time_vals);
        data{id} = data{id}(idx, :); % reorder all columns
        csvwrite(fullfile(datapath, 'video', filename{id}), data{id});
    end
end

% Get number of frames from the synchronized data

nframes{1} = size(data{1}, 1); nframes{2} = size(data{2}, 1);
fprintf('Face frames: %d\nEye frames: %d\n', nframes{1}, nframes{2});

%% Extract Tracking Coordinates
% Face cam
implant_x   = data{1}(:, 3:3:6);  implant_y   = data{1}(:, 4:3:7);
ear_x       = data{1}(:, 9:3:15); ear_y       = data{1}(:,10:3:16);
pupil_004_x = data{1}(:,18:3:39); pupil_004_y = data{1}(:,19:3:40);
eye_004_x   = data{1}(:,42:3:63); eye_004_y   = data{1}(:,43:3:64);
cheek_x     = data{1}(:,66:3:75); cheek_y     = data{1}(:,67:3:76);
nostril_x   = data{1}(:,78:3:87); nostril_y   = data{1}(:,79:3:88);
nose_x      = data{1}(:,[78,90:3:93]); nose_y  = data{1}(:,[79,91:3:94]);
jaw_x       = data{1}(:,96);      jaw_y       = data{1}(:,97);
tongue_x    = data{1}(:,99);      tongue_y    = data{1}(:,100);
spout_x    = data{1}(:,102);      spout_y    = data{1}(:,103); spout_likelihood = data{1}(:, 104);

% Eye cam
% general super model indexing
% pupil_007_x = data{2}(:,18:3:39); pupil_007_y = data{2}(:,19:3:40);
% eye_007_x   = data{2}(:,42:3:63); eye_007_y   = data{2}(:,43:3:64);

% small model indexing
pupil_007_x = data{2}(:,3:3:24); pupil_007_y = data{2}(:,4:3:25);
eye_007_x   = data{2}(:,27:3:48); eye_007_y   = data{2}(:,28:3:49);


center_body_parts = { ...
    'pupil_004', 'pupil_center_004'; ...
    'nostril',   'nostril_center'; ...
    'nose',      'nose_center'; ...
    'cheek',     'cheek_center'; ...
    'ear',       'ear_center'; ...
    };

%% Calculate Basic Metrics Per Frame
pupil_area_004 = nan(nframes{1},1); pupil_center_004 = nan(nframes{1},2);
eye_area_004 = nan(nframes{1},1);
nostril_area = nan(nframes{1},1); nostril_center = nan(nframes{1},2);
nose_area = nan(nframes{1},1); nose_center = nan(nframes{1},2);
cheek_center = nan(nframes{1},2); ear_center = nan(nframes{1},2);

for frame = 1:nframes{1}
    if any(~isnan(pupil_004_x(frame,:))) && any(~isnan(pupil_004_y(frame,:)))
        pupil_area_004(frame) = polyarea(pupil_004_x(frame,:), pupil_004_y(frame,:));
    end
    if any(~isnan(eye_004_x(frame,:))) && any(~isnan(eye_004_y(frame,:)))
        eye_area_004(frame) = polyarea(eye_004_x(frame,:), eye_004_y(frame,:));
    end
    if any(~isnan(nostril_x(frame,:))) && any(~isnan(nostril_y(frame,:)))
        nostril_area(frame) = polyarea(nostril_x(frame,:), nostril_y(frame,:));
    end
    if any(~isnan(nose_x(frame,:))) && any(~isnan(nose_y(frame,:)))
        nose_area(frame) = polyarea(nose_x(frame,:), nose_y(frame,:));
    end
    
    % Compute convex hull centers for each part
    for i = 1:size(center_body_parts,1)
        part = center_body_parts{i,1}; center_name = center_body_parts{i,2};
        x = eval([part '_x(frame,:)']); y = eval([part '_y(frame,:)']);
        valid_idx = ~isnan(x) & ~isnan(y);
        pts = [x(valid_idx)', y(valid_idx)'];
        if size(unique(pts,'rows'),1) < 3
            eval([center_name '(frame,:) = [NaN, NaN];']);
        else
            k = convhull(pts(:,1), pts(:,2));
            c = mean(pts(k,:),1);
            eval([center_name '(frame,:) = c;']);
        end
    end
end

jaw_center    = [jaw_x, jaw_y];
tongue_center = [tongue_x, tongue_y];
spout_center = [spout_x, spout_y];

%% Compute areas and center (eye cam)
pupil_area_007 = nan(nframes{2},1); pupil_center_007 = nan(nframes{2},2);
eye_area_007 = nan(nframes{2},1);
for frame = 1:nframes{2}
    if any(~isnan(pupil_007_x(frame,:))) && any(~isnan(pupil_007_y(frame,:)))
        pupil_area_007(frame) = polyarea(pupil_007_x(frame,:), pupil_007_y(frame,:));
        x = pupil_007_x(frame,:); y = pupil_007_y(frame,:);
        valid = ~isnan(x) & ~isnan(y);
        pts = [x(valid)', y(valid)'];
        if size(unique(pts,'rows'),1) < 3
            pupil_center_007(frame,:) = [NaN, NaN];
        else
            k = convhull(pts(:,1), pts(:,2));
            pupil_center_007(frame,:) = mean(pts(k,:),1);
        end
    end
    if any(~isnan(eye_007_x(frame,:))) && any(~isnan(eye_007_y(frame,:)))
        eye_area_007(frame) = polyarea(eye_007_x(frame,:), eye_007_y(frame,:));
    end
end

%% Compute Movement
marker_centers = {'pupil_center_004','nostril_center','nose_center','cheek_center','ear_center','jaw_center','tongue_center','spout_center'};
for i = 1:length(marker_centers)
    center = eval(marker_centers{i});
    mvt = [0; sqrt(sum(diff(center).^2,2))];
    eval([marker_centers{i} '_mvt = mvt;']);
end
% Pupil center from eye cam
pupil_center_007_mvt      = [0; sqrt(sum(diff(pupil_center_007).^2,2))];

%% Implant-centered normalization Face 004 markers normalisation on per-session implant width
implant_x_L = implant_x(:, 1); implant_x_R = implant_x(:, 2);
implant_y_L = implant_y(:, 1); implant_y_R = implant_y(:, 2);
midpoint = [(implant_x_L + implant_x_R)/2, (implant_y_L + implant_y_R)/2];
D_implant = sqrt((implant_x_R - implant_x_L).^2 + (implant_y_R - implant_y_L).^2);
D_implant_med  = nanmedian(D_implant);

for i = 1:length(marker_centers)
    C = eval(marker_centers{i});
    Cnorm = nan(size(C));
    for f = 1:nframes{1}
        if isnan(D_implant(f)) || D_implant(f)==0 || any(isnan(C(f,:))) || any(isnan(midpoint(f,:)))
            continue
        end
        Cnorm(f,:) = (C(f,:) - midpoint(f,:)) / D_implant(f);
    end
    eval([marker_centers{i} '_norm = Cnorm;']);
end

% ---- normalise areas, velocities, and movements by implant width ----------------
pupil_area_004_norm = pupil_area_004 ./ (D_implant_med.^2);   % ratio of areas
eye_area_004_norm = eye_area_004 ./ (D_implant_med.^2);
nose_area_norm = nose_area ./ (D_implant_med.^2);
nostril_area_norm = nostril_area ./ (D_implant_med.^2);

for i = 1:length(marker_centers)
    mvt = eval([marker_centers{i} '_mvt']);
    eval([marker_centers{i} '_mvt_norm = mvt ./ D_implant_med;']);
end

%% Eye 007 markers normalisation on per-session eye size
% Normalisation on the eye size
eye_x_L = eye_007_x(:, 1); eye_x_R = eye_007_x(:, 5);
eye_y_L = eye_007_y(:, 1); eye_y_R = eye_007_y(:, 5);

D_eye = sqrt((eye_x_R - eye_x_L).^2 + (eye_y_R - eye_y_L).^2);
D_eye_med  = nanmedian(D_eye);

pupil_area_007_norm = pupil_area_007 ./ (D_eye_med.^2);
eye_area_007_norm = eye_area_007 ./ (D_eye_med.^2);

pupil_center_007_norm = nan(size(pupil_center_007));

for f = 1:nframes{2}
    mid_eye = [(eye_x_L(f)+eye_x_R(f))/2, (eye_y_L(f)+eye_y_R(f))/2];
    if isnan(D_eye(f)) || D_eye(f)==0 || any(isnan(pupil_center_007(f,:)))
        continue
    end
    pupil_center_007_norm(f,:) = (pupil_center_007(f,:) - mid_eye) ./ D_eye(f);
end
pupil_center_007_mvt_norm = pupil_center_007_mvt ./ D_eye_med;

%% Convert to tsd
time_face = data{1}(:,1); time_eye = data{2}(:,1);
% Areas
pupil_area_004 = tsd(time_face, pupil_area_004_norm);
eye_area_004   = tsd(time_face, eye_area_004_norm);
nostril_area   = tsd(time_face, nostril_area_norm);
nose_area      = tsd(time_face, nose_area_norm);
spout_likelihood = tsd(time_face, spout_likelihood);

pupil_area_007 = tsd(time_eye, pupil_area_007_norm);
eye_area_007   = tsd(time_eye, eye_area_007_norm);

% Centers

for i = 1:length(marker_centers)
    eval([marker_centers{i} ' = tsd(time_face, ' marker_centers{i} '_norm);']);
    eval([marker_centers{i} '_mvt = tsd(time_face, ' marker_centers{i} '_mvt_norm);']);
end

pupil_center_007 = tsd(time_eye, pupil_center_007_norm);
pupil_center_007_mvt = tsd(time_eye, pupil_center_007_mvt_norm);

%% Save Processed Data
% Save the computed tsd objects to a MAT file in the DLC folder.

savePath = fullfile(datapath, 'video');

if ~exist(fullfile(savePath, 'DLC_data.mat'), 'file')
    save(fullfile(savePath, 'DLC_data.mat'), ...
        'Session_params','time_face','time_eye', ...
        'pupil_area_004','eye_area_004','nostril_area','nose_area', ...
        'pupil_area_007','eye_area_007', ...
        'pupil_center_004','pupil_center_007','nostril_center','nose_center','cheek_center','ear_center','jaw_center','tongue_center','spout_center',...
        'pupil_center_004_mvt','pupil_center_007_mvt','nostril_center_mvt','nose_center_mvt','cheek_center_mvt','ear_center_mvt','jaw_center_mvt','tongue_center_mvt','spout_center_mvt',...
        'spout_likelihood');
else
    save(fullfile(savePath, 'DLC_data.mat'), ...
        'Session_params','time_face','time_eye', ...
        'pupil_area_004','eye_area_004','nostril_area','nose_area', ...
        'pupil_area_007','eye_area_007', ...
        'pupil_center_004','pupil_center_007','nostril_center','nose_center','cheek_center','ear_center','jaw_center','tongue_center','spout_center',...
        'pupil_center_004_mvt','pupil_center_007_mvt','nostril_center_mvt','nose_center_mvt','cheek_center_mvt','ear_center_mvt','jaw_center_mvt','tongue_center_mvt','spout_center_mvt',...
        'spout_likelihood', '-append');
end

%% QC / Sanity check plots (time series + scatter) -------------------------
if isfield(Session_params,'do_plots') && Session_params.do_plots == 1
    
    % ---- time axis: auto-detect units (ts=1e-4s vs seconds)
    dt_face = nanmedian(diff(time_face));
    dt_eye  = nanmedian(diff(time_eye));
    
    if dt_face > 1
        t_face = time_face(:) / 1e4;   % seconds
    else
        t_face = time_face(:);         % already seconds
    end
    if dt_eye > 1
        t_eye = time_eye(:) / 1e4;
    else
        t_eye = time_eye(:);
    end
    
    qc_dir = fullfile(datapath, 'video');
    if ~exist(qc_dir, 'dir'), mkdir(qc_dir); end
    
    % helper indices for scatter
    ss_face = max(1, Session_params.qc_subsample);
    ss_eye  = max(1, Session_params.qc_subsample);
    idx_face = 1:ss_face:numel(t_face);
    idx_eye  = 1:ss_eye:numel(t_eye);
    
    %% Figure 1: Main dynamics (areas + movements + likelihood)
    f1 = figure('Visible', Session_params.fig_visibility);
    set(f1, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85]);
    
    try
        sgtitle([Session_params.animal_name '. Session: ' Session_params.session_selection], 'FontWeight', 'bold');
    catch
        suptitle([Session_params.animal_name '. Session: ' Session_params.session_selection]);
    end
    
    subplot(4,1,1);
    plot(t_face, pupil_area_004_norm); hold on;
    plot(t_face, eye_area_004_norm);
    plot(t_face, nostril_area_norm);
    plot(t_face, nose_area_norm);
    xlabel('Time (s)'); ylabel('Norm area');
    legend({'pupil004','eye004','nostril','nose'}, 'Location','best');
    title('Face cam: normalized areas'); grid on;
    
    subplot(4,1,2);
    plot(t_eye, pupil_area_007_norm); hold on;
    plot(t_eye, eye_area_007_norm);
    xlabel('Time (s)'); ylabel('Norm area');
    legend({'pupil007','eye007'}, 'Location','best');
    title('Eye cam: normalized areas'); grid on;
    
    subplot(4,1,3);
    plot(t_face, pupil_center_004_mvt_norm); hold on;
    plot(t_face, nostril_center_mvt_norm);
    plot(t_face, nose_center_mvt_norm);
    plot(t_face, jaw_center_mvt_norm);
    plot(t_face, tongue_center_mvt_norm);
    plot(t_face, spout_center_mvt_norm);
    xlabel('Time (s)'); ylabel('Norm movement');
    legend({'pupil','nostril','nose','jaw','tongue','spout'}, 'Location','best');
    title('Face cam: center movement (normalized)'); grid on;
    
    subplot(4,1,4);
    if exist('spout_likelihood','var') && ~isempty(spout_likelihood)
        sp = spout_likelihood;
        
        % ensure numeric vector
        if iscell(sp)
            try
                sp = cell2mat(sp);
            catch
                sp = [];
            end
        end
        if isstruct(sp) || isa(sp,'tsd')
            try
                sp = Data(sp);
            catch
                sp = [];
            end
        end
        sp = double(sp(:));
        
        % length match to face time
        n = min(numel(t_face), numel(sp));
        if n > 1
            plot(t_face(1:n), sp(1:n), 'k');
        end
    end
    xlabel('Time (s)'); ylabel('Likelihood');
    title('Spout likelihood (face cam)'); grid on;
    
    if Session_params.save_plots == 1
        print(f1, fullfile(qc_dir, 'QC_behaviour_timeseries.png'), '-dpng', '-r150');
    end
    
    
    %% Figure 3: Eye cam pupil center scatter + drift
    f3 = figure('Visible', Session_params.fig_visibility);
    set(f3, 'Units', 'Normalized', 'Position', [0.1 0.1 0.8 0.75]);
    
    subplot(2,1,1);
    scatter(pupil_center_007_norm(idx_eye,1), pupil_center_007_norm(idx_eye,2), 5, 'filled');
    axis equal; grid on;
    xlabel('x'); ylabel('y'); title('pupil center (eye cam)');
    
    subplot(2,1,2);
    plot(t_eye, pupil_center_007_norm(:,1)); hold on;
    plot(t_eye, pupil_center_007_norm(:,2));
    xlabel('Time (s)'); ylabel('Norm position');
    legend({'x','y'}, 'Location','best');
    title('pupil center drift (eye cam)'); grid on;
    
    if Session_params.save_plots == 1
        print(f3, fullfile(qc_dir, 'QC_behaviour_eye_pupil.png'), '-dpng', '-r150');
    end
    
    
    %% Figure 5: Mean bodypart coordinates (one plot)
    f5 = figure('Visible', Session_params.fig_visibility);
    set(f5, 'Units', 'Normalized', 'Position', [0.2 0.2 0.6 0.6]);
    hold on; grid on; axis equal;
    xlabel('x (norm)'); ylabel('y (norm)');
    title('Mean bodypart coordinates (centers)');
    
    % Face cam centers
    pts = {};
    if exist('pupil_center_004_norm','var') && ~isempty(pupil_center_004_norm)
        pts(end+1,:) = {'pupil_face', nanmean(pupil_center_004_norm(:,1)), nanmean(pupil_center_004_norm(:,2))};
    end
    if exist('nostril_center_norm','var') && ~isempty(nostril_center_norm)
        pts(end+1,:) = {'nostril', nanmean(nostril_center_norm(:,1)), nanmean(nostril_center_norm(:,2))};
    end
    if exist('nose_center_norm','var') && ~isempty(nose_center_norm)
        pts(end+1,:) = {'nose', nanmean(nose_center_norm(:,1)), nanmean(nose_center_norm(:,2))};
    end
    if exist('jaw_center_norm','var') && ~isempty(jaw_center_norm)
        pts(end+1,:) = {'jaw', nanmean(jaw_center_norm(:,1)), nanmean(jaw_center_norm(:,2))};
    end
    if exist('tongue_center_norm','var') && ~isempty(tongue_center_norm)
        pts(end+1,:) = {'tongue', nanmean(tongue_center_norm(:,1)), nanmean(tongue_center_norm(:,2))};
    end
    if exist('spout_center_norm','var') && ~isempty(spout_center_norm)
        pts(end+1,:) = {'spout', nanmean(spout_center_norm(:,1)), nanmean(spout_center_norm(:,2))};
    end
    
    % Eye cam pupil center
    if exist('pupil_center_007_norm','var') && ~isempty(pupil_center_007_norm)
        pts(end+1,:) = {'pupil_eye', nanmean(pupil_center_007_norm(:,1)), nanmean(pupil_center_007_norm(:,2))};
    end
    
    % Plot
    for i = 1:size(pts,1)
        name = pts{i,1};
        x = pts{i,2};
        y = pts{i,3};
        
        h = plot(x, y, 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
        text(x, y, ['  ' name], 'FontSize', 10, 'Color', get(h,'Color'));
    end
    
    legend({pts{:,1}}, 'Location','bestoutside');
    
    if Session_params.save_plots == 1
        print(f5, fullfile(qc_dir, 'QC_behaviour_mean_bodypart_coords.png'), '-dpng', '-r150');
    end
    
end






%% Legacy figure
% %% Plot figures
% f1 = figure('Visible', Session_params.fig_visibility);
% set(f1, 'Units', 'Normalized', 'Position', [0 0 1 1]);
%
% f2 = figure('Visible', Session_params.fig_visibility);
% set(f2, 'Units', 'Normalized', 'Position', [0 0 1 1]);
%
% f3 = figure('Visible', Session_params.fig_visibility);
% set(f3, 'Units', 'Normalized', 'Position', [0 0 1 1]);
%
% f4 = figure('Visible', Session_params.fig_visibility);
% set(f4, 'Units', 'Normalized', 'Position', [0 0 1 1]);
%
% f5 = figure('Visible', Session_params.fig_visibility);
% set(f5, 'Units', 'Normalized', 'Position', [0 0 1 1]);
%
% %% f1: Dynamics plots
% set(0, 'CurrentFigure', f1)
%
% try
%     sgtitle([Session_params.animal_name '. Session: '  Session_params.session_selection], 'FontWeight', 'bold')
% catch
%     suptitle([Session_params.animal_name '. Session: ' Session_params.session_selection], 'FontWeight', 'bold')
% end
% subplot(3,3,1)
% plot(data(:, 1), Data(pupil_004_x));
% title('pupil_x')
% xlabel('Time (s)')
%
% subplot(3,3,2)
% plot(data(:, 1), Data(pupil_004_y));
% title('pupil_y')
% xlabel('Time (s)')
%
% subplot(3,3,3)
% sparse_pupil_x = Data(pupil_004_x);
% sparse_pupil_y = Data(pupil_004_y);
% scatter(sparse_pupil_x(1:1000:end), sparse_pupil_y(1:1000:end));
% title('pupil')
% xlabel('Time (s)')
%
% % eye
% subplot(3,3,4)
% plot(data(:, 1), Data(eye_004_x))
% title('eye_x')
% xlabel('Time (s)')
%
% subplot(3,3,5)
% plot(data(:, 1), Data(eye_y))
% title('eye_y')
% xlabel('Time (s)')
%
% subplot(3,3,6)
% sparse_eye_x = Data(eye_004_x);
% sparse_eye_y = Data(eye_y);
% scatter(sparse_eye_x(1:1000:end), sparse_eye_y(1:1000:end))
% title('eye')
% xlabel('Time (s)')
%
% % another version of scatter that works in mobs
% % col_map=magma;
% %
% % figure
% % for i=1:8
% %     plot(eye_x(1:10:end,i), eye_y(1:10:end,i),'.','Color',col_map(32*i,:))
% %     hold on
% % end
% % for i=1:8
% %     plot(pupil_x(1:10:end,i), pupil_y(1:10:end,i),'.','Color',col_map(32*i,:))
% %     hold on
% % end
%
% % nostril
% subplot(3,3,7)
% plot(data(:, 1), Data(nostril_x))
% title('nostril_x')
% xlabel('Time (s)')
%
% subplot(3,3,8)
% plot(data(:, 1), Data(nostril_y))
% title('nostril_y')
% xlabel('Time (s)')
%
% subplot(3,3,9)
% sparse_nostril_x = Data(nostril_x);
% sparse_nostril_y = Data(nostril_y);
% scatter(sparse_nostril_x(1:1000:end), sparse_nostril_y(1:1000:end))
% title('nostril')
% xlabel('Time (s)')
%
% % another version of scatter that works in mobs
% % figure
% % for i=1:4
% %     plot(nostril_x(1:10:end,i), nostril_y(1:10:end,i),'.','Color',col_map(32*i,:))
% %     hold on
% % end

% %% f2: Areas plots
% set(0, 'CurrentFigure', f2)
% sgtitle([Session_params.animal_name '. Session: ' Session_params.session_selection], 'FontWeight', 'bold')
%
% subplot(311)
% plot(data(:, 1), Data(areas_pupil))
% hold on
% plot(data(:, 1), Data(areas_eye))
% title('pupil and eye areas evolution')
% xlabel('Time (s)')
%
% subplot(312)
% plot(data(:, 1), zscore(Data(areas_pupil)))
% hold on
% plot(data(:, 1), zscore(Data(areas_eye)))
% title('pupil and eye areas evolution (zscored)')
% xlabel('Time (s)')

% %% f3: Velocity/Acceleration plots
% set(0, 'CurrentFigure', f3)
% sgtitle([Session_params.animal_name '. Session: ' Session_params.session_selection], 'FontWeight', 'bold')
%
% subplot(411)
% plot(data(:, 1), Data(velocity_pupil_center));
% title('velocity of the pupil')
% legend({'velocity x', 'velocity y'})
% xlabel('Time (s)')
%
% subplot(412)
% plot(data(:, 1), Data(acceleration_pupil_center))
% title('Acceleration of the pupil')
% legend({'Acceleration x', 'Acceleration y'})
% xlabel('Time (s)')
%
% subplot(413)
% plot(data(:, 1), Data(velocity_nostril_center));
% title('velocity of the nostril')
% legend({'velocity x', 'velocity y'})
% xlabel('Time (s)')
%
% subplot(414)
% plot(data(:, 1), Data(acceleration_nostril_center))
% title('Acceleration of the nostril')
% legend({'Acceleration x', 'Acceleration y'})
% xlabel('Time (s)')

% %% f4: put everything together (eye + pupil)
% set(0, 'CurrentFigure', f4)
% sgtitle([Session_params.animal_name '. Session: ' Session_params.session_selection], 'FontWeight', 'bold')
%
% if exist('REMEpoch')
%     observed_rem_start = Start(REMEpoch);
%     observed_rem_end = End(REMEpoch);
% end
%
% % pupil center x
% axp1 = subplot(5, 1, 1);
% temp_data = Data(pupil_center);
% plot(data(:, 1), temp_data(:, 1))
% title('pupil center_x')
% xlabel('Time (s)')
%
% % pupil center y
% axp2 = subplot(5, 1, 2);
% plot(data(:, 1), temp_data(:, 2))
% title('pupil center_y')
% xlabel('Time (s)')
%
% % pupil and eye areas
% axp3 = subplot(5, 1, 3);
% plot(data(:, 1), zscore(Data(areas_pupil)))
% hold on
% plot(data(:, 1), zscore(Data(areas_eye)))
% title('pupil and eye areas (zscored)')
% legend({'pupil', 'eye'})
% xlabel('Time (s)')
%
% % pupil velocity x
% axp4 = subplot(5, 1, 4);
% temp_data = Data(velocity_pupil_center);
% plot(data(:, 1), temp_data(:, 1))
% title('pupil velocity_x')
% xlabel('Time (s)')
%
% % pupil velocity y
% axp5 = subplot(5, 1, 5);
% plot(data(:, 1), temp_data(:, 2))
% title('pupil velocity_y')
% xlabel('Time (s)')
%
% subplot_list = [axp1 axp2 axp3 axp4 axp5];
%
% % try
% %     load(fullfile(datapath, 'SleepScoring_OBGamma.mat'), 'REMEpoch')
% %     rem_start = Start(REMEpoch)/1e4;
% %     rem_end = End(REMEpoch)/1e4;
% % catch
% %     disp('No REM epoch found')
% % end
% % % plot REM episodes
% % if exist('REMEpoch')
% %     for i = 1:size(subplot_list, 2)
% %         axes(subplot_list(i))
% %         for j = 1:length(rem_start)
% %             patch([rem_start(j) rem_end(j) rem_end(j) rem_start(j)], [-5 -5 5 5], 'red', 'FaceAlpha', 0.2);
% %         end
% %         clear j
% %     end
% % end
%
% %% f5: put everything together (nostril)
% set(0, 'CurrentFigure', f5)
% sgtitle([Session_params.animal_name '. Session: ' Session_params.session_selection], 'FontWeight', 'bold')
%
% % nostril center x
% axn1 = subplot(5, 1, 1);
% temp_data = Data(nostril_center);
% plot(data(:, 1), temp_data(:, 1))
% title('nostril center_x')
%
% % nostril center y
% axn2 = subplot(5, 1, 2);
% plot(data(:, 1), temp_data(:, 2))
% title('nostril center_y')
%
% % nostril areas
% axn3 = subplot(5, 1, 3);
% plot(data(:, 1), Data(areas_nostril))
% title('nostril area')
%
% % nostril velocity x
% axn4 = subplot(5, 1, 4);
% temp_data = Data(velocity_nostril_center);
% plot(data(:, 1), temp_data(:, 1))
% title('nostril velocity_x')
%
% % nostril velocity y
% axn5 = subplot(5, 1, 5);
% plot(data(:, 1), temp_data(:, 2))
% title('nostril velocity_y')
%
% subplot_list = [axn1 axn2 axn3 axn4 axn5];
%
% % % plot REM episodes
% % if exist('REMEpoch')
% %     for i = 1:size(subplot_list, 2)
% %         axes(subplot_list(i))
% %         for j = 1:length(rem_start)
% %             patch([rem_start(j) rem_end(j) rem_end(j) rem_start(j)], [-5 -5 5 5], 'red', 'FaceAlpha', 0.2);
% %         end
% %         clear j
% %     end
% % end
%
% %% Save figures
% figFolder = fullfile(datapath, 'DLC/Figures');
% if ~exist(figFolder, 'dir')
%     mkdir(figFolder);
% end
%
% saveas(f1, fullfile(figFolder, [Session_params.animal_name '_eye_nostril_' Session_params.session_selection]), 'svg')
% saveas(f1, fullfile(figFolder, [Session_params.animal_name '_eye_nostril_' Session_params.session_selection]), 'png')
%
% saveas(f2, fullfile(figFolder, [Session_params.animal_name 'eye_nostril_areas_' Session_params.session_selection]), 'svg')
% saveas(f2, fullfile(figFolder, [Session_params.animal_name 'eye_nostril_areas_' Session_params.session_selection]), 'png')
%
% saveas(f3, fullfile(figFolder, [Session_params.animal_name 'eye_nostril_velocity_' Session_params.session_selection]), 'svg')
% saveas(f3, fullfile(figFolder, [Session_params.animal_name 'eye_nostril_velocity_' Session_params.session_selection]), 'png')
%
% saveas(f4, fullfile(figFolder, [Session_params.animal_name 'eye_pupil_all_' Session_params.session_selection]), 'svg')
% saveas(f4, fullfile(figFolder, [Session_params.animal_name 'eye_pupil_all_' Session_params.session_selection]), 'png')
%
% saveas(f5, fullfile(figFolder, [Session_params.animal_name 'nostril_all_' Session_params.session_selection]), 'svg')
% saveas(f5, fullfile(figFolder, [Session_params.animal_name 'nostril_all_' Session_params.session_selection]), 'png')

end