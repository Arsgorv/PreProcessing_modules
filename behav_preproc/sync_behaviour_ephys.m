function sync_behaviour_ephys(datapath, runTag)
%%%%%% Synchronization of video and OB (LFP) signal %%%%%%%%%
% This function:
%   1. Loads LFP data and extracts trigger timestamps.
%   2. Loads DLC tracking data.
%   3. Pairs DLC chunks to TTL trains conservatively.
%   4. Interpolates DLC to TTL count only when needed.
%   5. Saves synchronized csv + exact MAT sidecar times.
%
% Conservative updates:
%   - interpolation now uses index space (avoids duplicate x-values)
%   - exact synced times are saved in MAT sidecars
%   - RA experiment runs can be restricted to manifest phase windows
%   - catastrophic chunk mismatches are skipped cleanly with QC metadata

if nargin < 2, runTag = ''; end

Session_params.plt = 0;
gap_sec_video = 1.0;
max_allowed_mismatch = 200;
catastrophic_mismatch_abs = 2000;
catastrophic_mismatch_frac = 0.05;
TsRate = 1e4;

P = resolve_session_paths(datapath);
cfg = get_trigger_config(datapath);

if ~P.has_video
    error('sync_behaviour_ephys:NoVideoFolder', 'No video folder: %s', P.video_dir);
end
if ~P.has_lfp
    error('sync_behaviour_ephys:NoLFPData', 'No LFPData found for: %s', datapath);
end

outDir = P.video_dir;
if ~isempty(runTag)
    outDir = fullfile(P.video_dir, runTag);
    if ~exist(outDir,'dir'), mkdir(outDir); end
end

cam(1).name = 'FACE';
cam(1).role = 'face';
cam(1).ch   = cfg.face_ch;
cam(1).serial = '24934004';
cam(1).out_csv = fullfile(outDir, 'synchronized_DLC_data_face.csv');

cam(2).name = 'EYE';
cam(2).role = 'eye';
cam(2).ch   = cfg.eye_ch;
cam(2).serial = '24934007';
cam(2).out_csv = fullfile(outDir, 'synchronized_DLC_data_eye.csv');

dlc_all = dir(fullfile(P.video_dir, '*_filtered.csv'));
if isempty(dlc_all)
    error('sync_behaviour_ephys:NoDLC', 'No DLC csv found in %s', P.video_dir);
end

[runWin, hasRunWin] = get_run_window_if_needed(datapath, runTag);

for idx = 1:2
    if isempty(cam(idx).ch) || any(isnan(cam(idx).ch))
        warning('sync_behaviour_ephys:MissingChannel', ...
            '%s TTL channel undefined in get_trigger_config for %s', cam(idx).name, datapath);
        continue
    end

    lfp_file = fullfile(P.lfp_dir, ['LFP' num2str(cam(idx).ch) '.mat']);
    if ~exist(lfp_file, 'file')
        error('sync_behaviour_ephys:MissingLFPFile', 'Missing %s', lfp_file);
    end
    S = load(lfp_file);
    if ~isfield(S, 'LFP')
        error('sync_behaviour_ephys:BadLFPFile', 'Variable LFP not found in %s', lfp_file);
    end
    LFP_trig = S.LFP;

    time_lfp = Range(LFP_trig);
    LFP_trig_data = Data(LFP_trig);

    time_det = time_lfp(:);
    data_det = LFP_trig_data(:);

    if hasRunWin
        keepDet = time_det >= runWin.t0_ts & time_det <= runWin.t1_ts;
        if nnz(keepDet) > 10
            time_det = time_det(keepDet);
            data_det = data_det(keepDet);
        end
    end

    tmp = diff(data_det);
    if isempty(tmp) || all(tmp == 0)
        error('sync_behaviour_ephys:EmptyDiff', '%s: diff(LFP) empty/zero (ch %d)', cam(idx).name, cam(idx).ch);
    end

    threshold_diff = max(tmp) / 3;
    if threshold_diff > 0
        rise_fall_idx = find(tmp >= threshold_diff);
    else
        rise_fall_idx = [];
    end

    if numel(rise_fall_idx) < 2
        threshold_diff = max(abs(tmp)) / 3;
        rise_fall_idx = find(abs(tmp) >= threshold_diff);
    end

    threshold_event = 4;
    num_points = numel(rise_fall_idx);
    if num_points < 2
        error('sync_behaviour_ephys:TooFewEdges', '%s: too few TTL edges detected (ch %d).', cam(idx).name, cam(idx).ch);
    end

    event_indices = 1;
    for i = 2:num_points
        if (rise_fall_idx(i) - rise_fall_idx(i-1)) > threshold_event
            event_indices = [event_indices, i]; %#ok<AGROW>
        end
    end
    event_indices = [event_indices, num_points];

    num_events = numel(event_indices) - 1;
    peak_values  = nan(num_events, 1);
    peak_indices = nan(num_events, 1);

    for i = 1:num_events
        rr = rise_fall_idx(event_indices(i):event_indices(i+1)-1);
        event_data = LFP_trig_data(rr);
        if ~isempty(event_data)
            [~, max_idx] = max(event_data);
            peak_values(i)  = event_data(max_idx);
            peak_indices(i) = rr(max_idx);
        end
    end

    if any(isnan(peak_indices))
        idxNaN = find(isnan(peak_indices));
        peak_indices(idxNaN) = peak_indices(end-1) + round(nanmean(diff(peak_indices)));
        peak_values(idxNaN)  = peak_values(end-1)  + round(nanmean(diff(peak_values)));
        disp('NaNs in peak_indices/peak_values were replaced by approximations.');
    end

    interpeak_intervals = diff(peak_indices);
    if numel(interpeak_intervals) > 10
        z_scores = zscore(interpeak_intervals);
        threshold = 3.5;
        outlier_indices = find(z_scores > threshold | z_scores < -threshold);

        if ismember(1, outlier_indices)
            peak_indices(1) = [];
            peak_values(1)  = [];
            outlier_indices(outlier_indices == 1) = [];
        end

        if ~isempty(outlier_indices) && outlier_indices(end) >= numel(peak_indices)
            peak_indices(end) = [];
            peak_values(end)  = [];
        end
    end

    time_trig_all = time_det(peak_indices);
    fprintf('%s cam: first openephys trigger time: %fs (LFP time)\n', cam(idx).name, time_trig_all(1) / TsRate);

    time_trig_use = time_trig_all;
    usedPhaseWindow = hasRunWin;
    if isempty(time_trig_use)
        warning('sync_behaviour_ephys:EmptyPhaseWindow', ...
            '[%s] No TTLs remain after restricting to %s window.', cam(idx).name, runTag);
        continue
    end

    dlc_candidates = dlc_all(contains({dlc_all.name}, cam(idx).serial) & contains({dlc_all.name}, 'DLC'));
    if ~isempty(runTag)
        dlc_candidates = dlc_candidates(contains({dlc_candidates.name}, ['_' runTag]));
    end
    if isempty(dlc_candidates)
        warning('No DLC candidates for %s (%s) run=%s', cam(idx).name, cam(idx).serial, runTag);
        continue
    end

    if isempty(runTag)
        dlc_candidates = sort_training_dlc_candidates(dlc_candidates);
    else
        dlc_candidates = sort_by_suffix(dlc_candidates, runTag);
    end
    dlc_candidates = deduplicate_dlc_outputs(dlc_candidates);

    
    [~, ~, infoTr] = group_ttl_trains(time_trig_use, gap_sec_video);
    nTrains = infoTr.n_trials;
    trainIDs = 1:nTrains;

    if ~isempty(runTag) && ~usedPhaseWindow
        trainIDs = fallback_train_ids_for_run(dlc_all, dlc_candidates, cam(idx).serial, runTag, nTrains);
    end

    usedTrain = false(1, nTrains);
    allSync = [];
    allMeta = struct('dlc_name', {}, 'train_id', {}, 'nFrames', {}, 'nTTL', {}, ...
        'mismatch', {}, 'kept', {}, 'skip_reason', {}, 'runTag', {}, ...
        'usedPhaseWindow', {}, 'phase_t0_s', {}, 'phase_t1_s', {});

    for f = 1:numel(dlc_candidates)
        dlc_path = fullfile(P.video_dir, dlc_candidates(f).name);
        data = csvread(dlc_path, 3);
        nFrames = size(data, 1);

        bestScore = inf;
        bestTrain = NaN;
        bestTTL = NaN;

        for tr = trainIDs
            if tr < 1 || tr > nTrains, continue; end
            if usedTrain(tr), continue; end
            nTTL = sum(infoTr.trial_id == tr);
            score = abs(nTTL - nFrames);
            if score < bestScore
                bestScore = score;
                bestTrain = tr;
                bestTTL = nTTL;
            end
        end

        if ~isfinite(bestTrain)
            warning('[%s] No available TTL train for %s (run=%s)', cam(idx).name, dlc_candidates(f).name, runTag);
            continue
        end

        meta = struct();
        meta.dlc_name = dlc_candidates(f).name;
        meta.train_id = bestTrain;
        meta.nFrames = nFrames;
        meta.nTTL = bestTTL;
        meta.mismatch = bestScore;
        meta.kept = true;
        meta.skip_reason = '';
        meta.runTag = runTag;
        meta.usedPhaseWindow = usedPhaseWindow;
        if usedPhaseWindow
            meta.phase_t0_s = runWin.t0_ts / TsRate;
            meta.phase_t1_s = runWin.t1_ts / TsRate;
        else
            meta.phase_t0_s = NaN;
            meta.phase_t1_s = NaN;
        end

        skipThresh = max(catastrophic_mismatch_abs, round(catastrophic_mismatch_frac * max(nFrames, bestTTL)));
        if bestScore > skipThresh
            meta.kept = false;
            meta.skip_reason = sprintf('catastrophic mismatch (%d)', bestScore);
            warning('sync_behaviour_ephys:CatastrophicMismatch', ...
                '[%s %s] %s -> train %d mismatch=%d. Skipping this chunk.', ...
                cam(idx).name, runTag, dlc_candidates(f).name, bestTrain, bestScore);
            allMeta(end+1) = meta; %#ok<AGROW>
            continue
        elseif bestScore > max_allowed_mismatch
            warning('sync_behaviour_ephys:LargeMismatch', ...
                '[%s %s] %s -> train %d mismatch=%d. Keeping with warning.', ...
                cam(idx).name, runTag, dlc_candidates(f).name, bestTrain, bestScore);
        end

        usedTrain(bestTrain) = true;

        sel = (infoTr.trial_id == bestTrain);
        time_trig = time_trig_use(sel);

        [data_sync, savedIndices] = interpolate_dlc_to_ttl(data, time_trig);

        synchronizedData = [time_trig(:), data_sync];
        allSync = [allSync; synchronizedData]; %#ok<AGROW>
        allMeta(end+1) = meta; %#ok<AGROW>

        fprintf('[%s %s] %s -> train %d (TTL=%d frames=%d mismatch=%d)\n', ...
            cam(idx).name, runTag, dlc_candidates(f).name, bestTrain, bestTTL, nFrames, bestScore);

        if Session_params.plt == 1
            figure;
            tsec = time_trig(:) / TsRate;
            plot(tsec, data_sync(:,2), 'b-', 'LineWidth', 1.0); hold on;
            savedIndices = savedIndices(:);
            savedIndices = savedIndices(savedIndices >= 1 & savedIndices <= numel(tsec));
            plot(tsec(savedIndices), data_sync(savedIndices,2), 'ro', 'MarkerSize', 4);
            xlabel('Time (s)');
            ylabel('DLC dim #1');
            title([cam(idx).name ' DLC interpolation']);
            grid on;
        end
    end

    if isempty(allSync)
        warning('[%s] Empty synchronized output (run=%s)', cam(idx).name, runTag);
        if exist(cam(idx).out_csv, 'file')
            try, delete(cam(idx).out_csv); end %#ok<TRYNC>
        end
        write_sync_metadata_only(outDir, cam(idx).role, allMeta, runTag);
        continue
    end

    [~, o] = sort(allSync(:,1));
    allSync = allSync(o, :);

    csvwrite(cam(idx).out_csv, allSync);
    fprintf('Saved: %s\n', cam(idx).out_csv);

    allTime = allSync(:,1);
    sync_time_ts = allSync(:,1);
    syncMeta = allMeta;

    timesMat = strrep(cam(idx).out_csv, '.csv', '_times.mat');
    syncMat = strrep(cam(idx).out_csv, '.csv', '_sync.mat');
    legacyMat = fullfile(outDir, ['synchronized_DLC_times_' cam(idx).role '.mat']);

    save(timesMat, 'allTime', 'allMeta', 'runTag', '-v7.3');
    save(syncMat, 'sync_time_ts', 'syncMeta', 'runTag', '-v7.3');
    save(legacyMat, 'sync_time_ts', 'syncMeta', 'runTag', '-v7.3');

    make_sync_qc_figure(outDir, cam(idx), allMeta, runTag);

    fprintf('\n');
end
end

function [data_sync, savedIndices] = interpolate_dlc_to_ttl(data, time_trig)
nFrames = size(data, 1);
nTTL = numel(time_trig);
nDims = size(data, 2);

savedIndices = (1:nFrames)';
data_sync = data;

if nTTL == nFrames
    return
end

x_src = linspace(1, nTTL, nFrames);
x_q = 1:nTTL;

data_sync = nan(nTTL, nDims);
for col = 1:nDims
    data_sync(:, col) = interp1(x_src, data(:, col), x_q, 'linear', 'extrap');
end

savedIndices = round(x_src(:));
savedIndices(savedIndices < 1) = 1;
savedIndices(savedIndices > nTTL) = nTTL;
end

function trainIDs = fallback_train_ids_for_run(dlc_all, dlc_candidates, serial, runTag, nTrains)
trainIDs = 1:nTrains;

all_cam = dlc_all(contains({dlc_all.name}, serial) & contains({dlc_all.name}, 'DLC'));
nCond = sum(contains({all_cam.name}, '_Conditioning'));
nPost = sum(contains({all_cam.name}, '_PostTest'));

if strcmpi(runTag, 'Conditioning')
    trainIDs = 1:min(nCond, nTrains);
elseif strcmpi(runTag, 'PostTest')
    a = min(nCond, nTrains) + 1;
    b = min(nCond + nPost, nTrains);
    trainIDs = a:b;
end

if isempty(trainIDs)
    trainIDs = 1:nTrains;
end
end

function d = sort_by_suffix(d, runTag)
if isempty(d), return; end
k = zeros(numel(d),1);
for i = 1:numel(d)
    nm = d(i).name;
    tok = regexp(nm, ['_' runTag '_(\d+)'], 'tokens', 'once');
    if isempty(tok)
        k(i) = 0;
    else
        k(i) = str2double(tok{1});
    end
end
[~, ord] = sort(k);
d = d(ord);
end

function d = sort_training_dlc_candidates(d)
if isempty(d), return; end
nm = {d.name};
[~, ord] = sort(lower(nm));
d = d(ord);
end

function d = deduplicate_dlc_outputs(d)
if isempty(d), return; end
keep = true(numel(d),1);
keys = cell(numel(d),1);
for i = 1:numel(d)
    nm = d(i).name;
    tok = regexp(nm, '^(.*)DLC_', 'tokens', 'once');
    if isempty(tok)
        keys{i} = lower(nm);
    else
        keys{i} = lower(tok{1});
    end
end
u = unique(keys, 'stable');
out = d([]);
for i = 1:numel(u)
    idx = find(strcmp(keys, u{i}));
    if numel(idx) == 1
        out(end+1) = d(idx); %#ok<AGROW>
    else
        [~, ord] = sort([d(idx).datenum], 'descend');
        out(end+1) = d(idx(ord(1))); %#ok<AGROW>
    end
end
d = out;
end

function [runWin, hasRunWin] = get_run_window_if_needed(datapath, runTag)
runWin = struct('t0_ts', NaN, 't1_ts', NaN);
hasRunWin = false;

if isempty(runTag)
    return
end

if ~strcmpi(runTag, 'Conditioning') && ~strcmpi(runTag, 'PostTest')
    return
end

mf = fullfile(datapath, 'analysis', 'run_manifest_RAexp.mat');
try
    if exist(mf, 'file')
        S = load(mf, 'RUN');
        RUN = S.RUN;
    else
        RUN = RAE_make_run_manifest(datapath);
    end

    if isfield(RUN, 'run') && isfield(RUN.run, runTag) && ...
            isfield(RUN.run.(runTag), 't0_ts') && isfield(RUN.run.(runTag), 't1_ts') && ...
            isfinite(RUN.run.(runTag).t0_ts) && isfinite(RUN.run.(runTag).t1_ts)
        runWin.t0_ts = RUN.run.(runTag).t0_ts;
        runWin.t1_ts = RUN.run.(runTag).t1_ts;
        hasRunWin = true;
    end
catch ME
    warning('sync_behaviour_ephys:RunManifest', ...
        'Could not load/create run manifest for %s (%s). Falling back to legacy train allocation.', ...
        runTag, ME.message);
end
end

function write_sync_metadata_only(outDir, camRole, allMeta, runTag)
if nargin < 3 || isempty(allMeta)
    return
end
metaMat = fullfile(outDir, ['synchronized_DLC_times_' camRole '.mat']);
syncMeta = allMeta; %#ok<NASGU>
sync_time_ts = []; %#ok<NASGU>
save(metaMat, 'sync_time_ts', 'syncMeta', 'runTag', '-v7.3');
end

function make_sync_qc_figure(outDir, cam, allMeta, runTag)
try
    if isempty(allMeta)
        return
    end

    n = numel(allMeta);
    nFrames = nan(n,1);
    nTTL = nan(n,1);
    mismatch = nan(n,1);
    kept = false(n,1);

    for i = 1:n
        if isfield(allMeta(i), 'nFrames'), nFrames(i) = allMeta(i).nFrames; end
        if isfield(allMeta(i), 'nTTL'), nTTL(i) = allMeta(i).nTTL; end
        if isfield(allMeta(i), 'mismatch'), mismatch(i) = allMeta(i).mismatch; end
        if isfield(allMeta(i), 'kept'), kept(i) = allMeta(i).kept; end
    end

    f = figure('Visible', 'off', 'Color', 'w');
    set(f, 'Units', 'Pixels', 'Position', [100 100 1200 500]);

    subplot(1,2,1);
    bar([nFrames nTTL], 'grouped');
    xlabel('Chunk');
    ylabel('Count');
    title([cam.name ' sync counts']);
    legend({'DLC frames','TTL pulses'}, 'Location', 'best');
    set(gca, 'XTick', 1:n);
    grid on;

    subplot(1,2,2);
    stem(1:n, mismatch, 'filled');
    hold on;
    bad = find(~kept);
    if ~isempty(bad)
        plot(bad, mismatch(bad), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    end
    xlabel('Chunk');
    ylabel('|TTL - frames|');
    title([cam.name ' sync mismatch']);
    set(gca, 'XTick', 1:n);
    grid on;

    tag = '';
    if ~isempty(runTag)
        tag = ['_' runTag];
    end
    saveas(f, fullfile(outDir, ['QC_sync_' lower(cam.role) tag '.png']));
    close(f);
catch ME
    warning('sync_behaviour_ephys:QCPlotFailed', ...
        'Could not save sync QC figure for %s (%s)', cam.name, ME.message);
end
end
