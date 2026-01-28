function QC = RA_sanity_check_datasync_active(datapath, varargin)
% RA_sanity_check_datasync_active
% Sanity figures for Master_data_sync_preproc (RA).
% Produces 2 figures per session:
%   Fig1: global alignment overview
%   Fig2: behavior + spout QC
%
% QC output contains detection summary for spout.

% -------------------- defaults --------------------
SaveFig = true;
OutDir = fullfile(datapath, 'stim', 'sanity');
FigVisible = 'on';

ArrLagNom = 3.95;     % seconds, expected spout arrival (trial-relative)
ArrLagTol = 0.50;     % seconds
WinPad = 1.5;         % seconds for spout QC window around ArrLagNom

SpoutThrHigh = 0.60;  % "good" threshold
SpoutThrLow  = 0.30;  % fallback threshold for max-based pick
SmoothN = 20;         % samples for runmean on likelihood

MaxTrialsHeatmap = 120;  % limit for likelihood heatmap
NExampleTraces = 12;     % example trials shown as traces

% -------------------- parse varargin --------------------
if ~isempty(varargin)
    for k = 1:2:numel(varargin)
        key = varargin{k};
        val = varargin{k+1};
        switch lower(key)
            case 'savefig'
                SaveFig = logical(val);
            case 'outdir'
                OutDir = val;
            case 'figvisible'
                FigVisible = val;

            case 'arrlagnom'
                ArrLagNom = val;
            case 'arrlagtol'
                ArrLagTol = val;
            case 'winpad'
                WinPad = val;

            case 'spoutthrhigh'
                SpoutThrHigh = val;
            case 'spoutthrlow'
                SpoutThrLow = val;
            case 'smoothn'
                SmoothN = val;

            case 'maxtrialsheatmap'
                MaxTrialsHeatmap = val;
            case 'nexampletraces'
                NExampleTraces = val;

            otherwise
                error('Unknown parameter: %s', key);
        end
    end
end

% -------------------- load Master_sync + Baphy --------------------
msfile = fullfile(datapath, 'Master_sync.mat');
if ~exist(msfile,'file')
    error('Master_sync.mat not found: %s', msfile);
end
S = load(msfile);

if isfield(S,'Baphy')
    B = S.Baphy;
else
    bfile = fullfile(datapath, 'stim', 'Baphy_RA.mat');
    if exist(bfile,'file')
        Sb = load(bfile);
        if isfield(Sb,'Baphy')
            B = Sb.Baphy;
        else
            error('Baphy not found in %s', bfile);
        end
    else
        error('Baphy not found in Master_sync.mat and %s missing', bfile);
    end
end

if ~isfield(S,'Epochs')
    if isfield(S,'trigOE')
        Epochs = RA_build_epochs_active(datapath, S.trigOE, B);
    else
        error('Epochs missing in Master_sync.mat and trigOE missing to rebuild.');
    end
else
    Epochs = S.Epochs;
end

if isfield(S,'trigOE')
    trigOE = S.trigOE;
else
    trigOE = struct();
end

if ~exist(OutDir,'dir')
    mkdir(OutDir);
end

[~,sessname] = fileparts(datapath);

% -------------------- assemble key times --------------------
TsRate = 1e4;

t_trial0 = B.trial.abs_trialstart_s(:);
t_trial1 = B.trial.abs_trialend_s(:);

t_stim0  = B.trial.abs_stim_start_s(:);
t_stim1  = B.trial.abs_stim_stop_s(:);

t_arr    = B.trial.spout_arrival_abs_s(:);

n = B.n_trials;

isRef = strcmpi(B.trial.type(1:n), 'Reference');
isTar = strcmpi(B.trial.type(1:n), 'Target');

isNoSound = false(n,1);
if isfield(B.trial,'is_nosound'), isNoSound = logical(B.trial.is_nosound(1:n)); end

isNoMotorExclusive = false(n,1);
if isfield(B.trial,'is_nomotor')
    isNoMotorExclusive = logical(B.trial.is_nomotor(1:n));
end

isNoMotorRaw = isNoMotorExclusive;
if isfield(B,'idx') && isfield(B.idx,'NoMotor') && ~isempty(B.idx.NoMotor)
    isNoMotorRaw = false(n,1);
    nm = B.idx.NoMotor(:);
    nm = nm(nm>=1 & nm<=n);
    isNoMotorRaw(nm) = true;
end

% -------------------- FIGURE 1: Global alignment overview --------------------
f1 = figure('Color','w', 'Visible', FigVisible);
set(f1, 'Units','pixels', 'Position',[50 50 1800 900]);

ax1 = gobjects(6,1);

% Panel 1: Event raster (trial/stim + masks + arrival) [minutes]
ax1(1) = subplot(3,2,1); hold on;
x0 = t_trial0/60; x1 = t_trial1/60;
xs = t_stim0/60;  xe = t_stim1/60;
xa = t_arr/60;

y_trial = 9; y_stim = 8;
y_tar = 7; y_ref = 6;
y_ns  = 5; y_nm = 4; y_nmx = 3;
y_arr = 2;

m = isfinite(t_trial0) & isfinite(t_trial1) & (t_trial1>t_trial0);
plot([x0(m) x1(m)]', [y_trial y_trial]'*ones(1,sum(m)), '-', 'LineWidth',1);

m = isfinite(t_stim0) & isfinite(t_stim1) & (t_stim1>t_stim0);
plot([xs(m) xe(m)]', [y_stim y_stim]'*ones(1,sum(m)), '-', 'LineWidth',1);

m = isfinite(t_stim0) & isRef;
plot(xs(m), y_ref*ones(sum(m),1), '.', 'MarkerSize',10);
m = isfinite(t_stim0) & isTar;
plot(xs(m), y_tar*ones(sum(m),1), '.', 'MarkerSize',10);

m = isNoSound & isfinite(t_trial0);
plot(x0(m), y_ns*ones(sum(m),1), 'x', 'MarkerSize',6);
m = isNoMotorRaw & isfinite(t_trial0);
plot(x0(m), y_nm*ones(sum(m),1), '*', 'MarkerSize',6);
m = isNoMotorExclusive & isfinite(t_trial0);
plot(x0(m), y_nmx*ones(sum(m),1), 'o', 'MarkerSize',4);

m = isfinite(t_arr);
plot(xa(m), y_arr*ones(sum(m),1), '.', 'MarkerSize',10);

yticks([y_arr y_nmx y_nm y_ns y_ref y_tar y_stim y_trial]);
yticklabels({'arrival','nomotorX','nomotor','nosound','ref stim','tar stim','stim','trial'});
xlabel('Time (min)'); title('Event raster');
grid on; box on;

% Panel 2: Trial/stim duration diagnostics [trial index]
ax1(2) = subplot(3,2,2); hold on;
trial_dur = t_trial1 - t_trial0;
stim_dur  = t_stim1 - t_stim0;

plot(trial_dur, '.', 'MarkerSize',8);
plot(stim_dur,  '.', 'MarkerSize',8);
ylabel('Duration (s)'); xlabel('Trial index');
legend({'trial dur','stim dur'}, 'Location','best');
title('Durations vs trial index');
grid on; box on;

% Panel 3: Stim position inside trial [trial index]
ax1(3) = subplot(3,2,3); hold on;
stim_on_rel  = B.trial.rel_stim_start(1:n);
stim_off_rel = B.trial.rel_stim_stop(1:n);

plot(stim_on_rel,  '.', 'MarkerSize',8);
plot(stim_off_rel, '.', 'MarkerSize',8);
ylabel('Time from trial start (s)'); xlabel('Trial index');
legend({'stim on','stim off'}, 'Location','best');
title('Stim timing (trial-relative)');
grid on; box on;

% Panel 4: Overlay representative streams with stim markers [minutes]
ax1(4) = subplot(3,2,4); hold on;
plotted_any = false;

if isfield(S,'cat_tsd') && isfield(S.cat_tsd,'data')
    try
        fus = S.cat_tsd.data;
        fus_t = Range(fus,'s')/60;
        fus_y = Data(fus);
        fus_y = reshape(fus_y', S.cat_tsd.Nx, S.cat_tsd.Ny, size(fus_y,1));
        fus_y = squeeze(mean(fus_y, [1 2]));
        plot(fus_t, fus_y, '-', 'LineWidth',1);
        plotted_any = true;
    catch
    end
else
    fusdir = fullfile(datapath,'fUS');
    if exist(fusdir,'dir')
        d = dir(fullfile(fusdir,'*slice_A.mat'));
        if ~isempty(d)
            try
                F = load(fullfile(fusdir, d(1).name));
                if isfield(F,'cat_tsd') && isfield(F.cat_tsd,'data')
                    fus = F.cat_tsd.data;
                    fus_t = Range(fus,'s')/60;
                    fus_y = Data(fus);
                    fus_y = reshape(fus_y', F.cat_tsd.Nx, F.cat_tsd.Ny, size(fus_y,1));
                    fus_y = squeeze(mean(fus_y, [1 2]));
                    plot(fus_t, fus_y, '-', 'LineWidth',1);
                    plotted_any = true;
                end
            catch
            end
        end
    end
end

if isfield(S,'LFP')
    try
        lfp = S.LFP;
        lfp_t = Range(lfp,'s')/60;
        lfp_y = Data(lfp);
        plot(lfp_t, lfp_y, '-', 'LineWidth',1);
        plotted_any = true;
    catch
    end
else
    lfppath = fullfile(datapath,'ephys','LFPData');
    if exist(lfppath,'dir')
        d = dir(fullfile(lfppath,'LFP*.mat'));
        if ~isempty(d)
            try
                L = load(fullfile(lfppath, 'LFP22'));
                if isfield(L,'LFP')
                    lfp = L.LFP;
                    lfp_t = Range(lfp,'s')/60;
                    lfp_y = Data(lfp);
                    plot(lfp_t, lfp_y, '-', 'LineWidth',1);
                    plotted_any = true;
                end
            catch
            end
        end
    end
end

m = isfinite(t_stim0);
yl = ylim;
for i = find(m(:))'
    x = t_stim0(i)/60;
    line([x x], yl, 'LineStyle',':');
end
ylim(yl);

xlabel('Time (min)'); title('Representative streams + stim onsets');
grid on; box on;
if ~plotted_any
    text(0.1,0.5,'Could not load fUS/LFP streams here.','Units','normalized');
end

% Panel 5: OE trial starts vs Baphy abs trial starts [trial index]
ax1(5) = subplot(3,2,5); hold on;
if isfield(trigOE,'baphy') && isfield(trigOE.baphy,'trial_start_ts') && ~isempty(trigOE.baphy.trial_start_ts)
    oe_ts = trigOE.baphy.trial_start_ts(:);
    m2 = min(numel(oe_ts), n);
    b_ts = round(B.trial.abs_trialstart_s(1:m2) * TsRate);
    d = double(b_ts) - double(oe_ts(1:m2));
    plot(d, '.', 'MarkerSize',8);
    ylabel('Baphy_ts - OE_ts (ticks)'); xlabel('Trial index');
    title('Start alignment check (ticks)');
    grid on; box on;
else
    axis off;
    text(0.1,0.5,'No trigOE.baphy.trial_start_ts to compare.','Units','normalized');
end

% Panel 6: Coverage bars per stream [minutes]
ax1(6) = subplot(3,2,6); hold on;
labels = {};
y = 1;

labels{end+1} = 'Baphy(trials)';
plot([min(t_trial0)/60 max(t_trial1)/60], [y y], '-', 'LineWidth',6); y=y+1;

if exist('lfp','var')
    labels{end+1} = 'LFP';
    plot([min(Range(lfp,'s'))/60 max(Range(lfp,'s'))/60], [y y], '-', 'LineWidth',6); y=y+1;
end
if exist('fus','var')
    labels{end+1} = 'fUS';
    plot([min(Range(fus,'s'))/60 max(Range(fus,'s'))/60], [y y], '-', 'LineWidth',6); y=y+1;
end

if isfield(trigOE,'face_cam') && isfield(trigOE.face_cam,'time_s') && ~isempty(trigOE.face_cam.time_s)
    labels{end+1} = 'face_cam';
    plot([min(trigOE.face_cam.time_s)/60 max(trigOE.face_cam.time_s)/60], [y y], '-', 'LineWidth',6); y=y+1;
end
if isfield(trigOE,'eye_cam') && isfield(trigOE.eye_cam,'time_s') && ~isempty(trigOE.eye_cam.time_s)
    labels{end+1} = 'eye_cam';
    plot([min(trigOE.eye_cam.time_s)/60 max(trigOE.eye_cam.time_s)/60], [y y], '-', 'LineWidth',6); y=y+1;
end

yticks(1:(y-1));
yticklabels(labels);
xlabel('Time (min)'); title('Stream coverage');
grid on; box on;

sgtitle(sprintf('RA sync sanity — %s', strrep(sessname,'_','\_')));

% ---- link x-axes (requested) ----
linkaxes(ax1([1 4 6]), 'x');   % minutes panels
linkaxes(ax1([2 3 5]), 'x');   % trial-index panels

if SaveFig
    saveas(f1, fullfile(OutDir, sprintf('%s_sanity1_global.png', sessname)));
end

% -------------------- FIGURE 2: Behavior + spout QC --------------------
f2 = figure('Color','w', 'Visible', FigVisible);
set(f2, 'Units','pixels', 'Position',[70 70 1800 900]);

QC = struct();
QC.datapath = datapath;
QC.sessname = sessname;
QC.n_trials = n;

arr_rel = B.trial.spout_arrival_rel_s(1:n);

% Panel 1: Spout arrival relative timing (trial-relative) + ref/tar
subplot(2,3,1); hold on;
histogram(arr_rel(isRef & isfinite(arr_rel)), 30);
histogram(arr_rel(isTar & isfinite(arr_rel)), 30);
xlabel('Arrival - trial start (s)'); ylabel('count');
legend({'Ref','Tar'}, 'Location','best');
title('Spout arrival (trial-relative)');
grid on; box on;

% Panel 2: [stim onset/offset -> arrival] duration distributions
subplot(2,3,2); hold on;
d_on  = B.trial.spout_arrival_abs_s(1:n) - B.trial.abs_stim_start_s(1:n);
d_off = B.trial.spout_arrival_abs_s(1:n) - B.trial.abs_stim_stop_s(1:n);
histogram(d_on(isfinite(d_on)), 40);
histogram(d_off(isfinite(d_off)), 40);
xlabel('Duration (s)'); ylabel('count');
legend({'stim on -> arrival','stim off -> arrival'}, 'Location','best');
title('Stim-to-arrival durations');
grid on; box on;

% Panel 3: Motor vs arrival (fixed masking)
subplot(2,3,3); hold on;
if isfield(Epochs,'motor_move_start_ts') && isfield(Epochs,'motor_arrival_ts') ...
        && ~isempty(Epochs.motor_move_start_ts) && ~isempty(Epochs.motor_arrival_ts)

    mot0 = double(Epochs.motor_move_start_ts(:))/TsRate;
    tr0  = B.trial.abs_trialstart_s(1:n);

    mot0_rel = mot0 - tr0;

    m = isfinite(mot0_rel) & isfinite(arr_rel);
    plot(mot0_rel(m), arr_rel(m), '.', 'MarkerSize',10);

    xlabel('motor start rel (s)');
    ylabel('spout arrival rel (s)');
    title('Motor start vs spout arrival');
    grid on; box on;
else
    axis off;
    text(0.1,0.5,'Motor times not available in Epochs.','Units','normalized');
end

% Panel 4: Spout likelihood heatmap aligned to trial start (interp1 unique fix)
subplot(2,3,4); hold on;

dlcfile = fullfile(datapath, 'video', 'DLC_data.mat');
if exist(dlcfile,'file')
    D = load(dlcfile, 'spout_likelihood');
    if isfield(D,'spout_likelihood') && ~isempty(D.spout_likelihood)
        spL = D.spout_likelihood;

        sp = [];
        try
            temp = runmean(Data(spL), SmoothN);
            sp = tsd(Range(spL), temp);
        catch
        end

        if ~isempty(sp)
            t_s = Range(sp,'s');
            y_s = Data(sp);

            % ensure unique time stamps globally
            [t_su, iau] = unique(t_s, 'stable');
            y_su = y_s(iau);
            t_s = t_su;
            y_s = y_su;

            winBeg = ArrLagNom - WinPad;
            winEnd = ArrLagNom + WinPad;

            okTr = find(isfinite(B.trial.abs_trialstart_s(1:n)));
            if numel(okTr) > MaxTrialsHeatmap
                okTr = okTr(round(linspace(1,numel(okTr),MaxTrialsHeatmap)));
            end

            dt = median(diff(t_s));
            if ~isfinite(dt) || dt<=0
                dt = 1/200;
            end
            tt = (winBeg:dt:winEnd)';
            L = nan(numel(okTr), numel(tt));

            maxL = nan(numel(okTr),1);
            tMax = nan(numel(okTr),1);

            for i = 1:numel(okTr)
                tr = okTr(i);
                t0 = B.trial.abs_trialstart_s(tr);
                a = t0 + winBeg;
                b = t0 + winEnd;

                idx = find(t_s >= a & t_s <= b);
                if numel(idx) < 5, continue; end

                tv = t_s(idx) - t0;
                lv = y_s(idx);

                % unique within-trial window for interp1
                [tvu, ia] = unique(tv, 'stable');
                lvu = lv(ia);

                if numel(tvu) < 5, continue; end

                L(i,:) = interp1(tvu, lvu, tt, 'linear', NaN);

                [mx, imx] = max(lvu);
                maxL(i) = mx;
                tMax(i) = tvu(imx);
            end

            imagesc(tt, 1:size(L,1), L);
            set(gca,'YDir','normal');
            xlabel('time from trial start (s)');
            ylabel('trial (subset)');
            title('Spout likelihood aligned (window)');
            colorbar;
            grid on; box on;

            detHigh = maxL >= SpoutThrHigh;
            detLow  = maxL >= SpoutThrLow;

            QC.spout.maxL = maxL;
            QC.spout.tMax_rel = tMax;
            QC.spout.detHigh_frac = mean(detHigh(isfinite(maxL)));
            QC.spout.detLow_frac  = mean(detLow(isfinite(maxL)));
            QC.spout.n_heatmap = numel(okTr);

        else
            axis off;
            text(0.1,0.5,'Could not build tsd for spout_likelihood.','Units','normalized');
        end
    else
        axis off;
        text(0.1,0.5,'DLC_data.mat present but spout_likelihood missing/empty.','Units','normalized');
    end
else
    axis off;
    text(0.1,0.5,'DLC_data.mat not found.','Units','normalized');
end

% Panel 5: Example likelihood traces + Baphy arrival picks
subplot(2,3,5); hold on;
if exist('sp','var') && ~isempty(sp)
    okTr = find(isfinite(B.trial.abs_trialstart_s(1:n)));
    if ~isempty(okTr)
        if numel(okTr) > NExampleTraces
            rp = randperm(numel(okTr), NExampleTraces);
            okTr = okTr(rp);
        end

        t_s = Range(sp,'s');
        y_s = Data(sp);

        [t_su, iau] = unique(t_s, 'stable');
        y_su = y_s(iau);
        t_s = t_su;
        y_s = y_su;

        winBeg = ArrLagNom - WinPad;
        winEnd = ArrLagNom + WinPad;

        for i = 1:numel(okTr)
            tr = okTr(i);
            t0 = B.trial.abs_trialstart_s(tr);
            a = t0 + winBeg;
            b = t0 + winEnd;

            idx = find(t_s >= a & t_s <= b);
            if numel(idx) < 5, continue; end

            tv = t_s(idx) - t0;
            lv = y_s(idx);

            [tvu, ia] = unique(tv, 'stable');
            lvu = lv(ia);

            plot(tvu, lvu);

            ar = B.trial.spout_arrival_rel_s(tr);
            if isfinite(ar)
                plot(ar, 0.95, 'kv', 'MarkerFaceColor','k', 'MarkerSize',4);
            end
        end

        yl = ylim;
        line([ArrLagNom-ArrLagTol ArrLagNom-ArrLagTol], yl, 'LineStyle',':');
        line([ArrLagNom+ArrLagTol ArrLagNom+ArrLagTol], yl, 'LineStyle',':');
        ylim(yl);

        xlabel('time from trial start (s)');
        ylabel('spout likelihood');
        title('Example likelihood traces (window)');
        grid on; box on;
    end
else
    axis off;
    text(0.1,0.5,'No spout likelihood available for traces.','Units','normalized');
end

% Panel 6: Hit rate / performance vs time with masks
subplot(2,3,6); hold on;
if isfield(B,'Perf') && isfield(B.Perf,'hit_rate_all') && ~isempty(B.Perf.hit_rate_all)
    hr = B.Perf.hit_rate_all(1:n);
    plot(t_trial0/60, hr, '.', 'MarkerSize',10);

    m = isNoSound & isfinite(hr);
    plot(t_trial0(m)/60, hr(m), 'x', 'MarkerSize',8);

    m = isNoMotorExclusive & isfinite(hr);
    plot(t_trial0(m)/60, hr(m), 'o', 'MarkerSize',4);

    xlabel('Time (min)'); ylabel('Hit rate');
    legend({'all','nosound','nomotorX'}, 'Location','best');
    title('Performance with masks');
    grid on; box on;
else
    axis off;
    text(0.1,0.5,'B.Perf.hit_rate_all not available.','Units','normalized');
end

sgtitle(sprintf('RA behavior + spout QC — %s', strrep(sessname,'_','\_')));

if SaveFig
    saveas(f2, fullfile(OutDir, sprintf('%s_sanity2_behavior_spout.png', sessname)));
end

% -------------------- spout summary --------------------
QC.spout.n_trials_total = n;
QC.spout.n_arrival_finite = sum(isfinite(B.trial.spout_arrival_abs_s(1:n)));
QC.spout.frac_arrival_finite = QC.spout.n_arrival_finite / n;

QC.spout.recommendation = { ...
    'Keep per-trial spout arrival as NaN when likelihood evidence is insufficient (do not interpolate as if it were continuous).', ...
    'For analyses needing a single timing per condition: use median over detected trials (Ref/Tar separately) and keep a mask of detected vs missing.', ...
    'If you must impute per trial: fill missing with condition median ONLY for visualization/alignment and keep an "imputed" flag to exclude from stats.'};

end