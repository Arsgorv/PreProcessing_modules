function QC = FMA_epoch_sanity_figures(datapath, wcRoot, RUN, opts)
% FMA_epoch_sanity_figures
% Multi-column figure: one column per high-level epoch.
%
% Rows:
% 1) Trial-aligned raster (representative unit; head-fixed only)
% 2) Trial-aligned PSTH (representative unit; head-fixed only)
% 3) Population raster across time (relative to epoch start)
% 4) Binned raster heatmap (MakeQfromS)
% 5) "STRF-like": SoundName x time PSTH map (rep unit; head-fixed only)
% 6) ZETA summary (head-fixed only)
%
% Saves:
%   <wcRoot>/sanity/FMA_epoch_sanity.png (+ .mat with QC)

if nargin < 4, opts = struct(); end
if ~isfield(opts,'phases'), opts.phases = {'PreSleep','Conditioning','PostSleep','PostTest'}; end

TsRate = 1e4;

outDir = fullfile(wcRoot, 'sanity');
if ~exist(outDir,'dir'), mkdir(outDir); end

nCols = numel(opts.phases);
nRows = 6;

fh = figure('Color','w','Position',[100 80 1850 1000]);

QC = struct();
QC.phases = opts.phases;

for c = 1:nCols
    phase = opts.phases{c};
    phDir = fullfile(wcRoot, phase);
    uFile = fullfile(phDir, ['FMA_units_' phase '.mat']);

    % epoch bounds
    t0_ts = NaN; t1_ts = NaN;
    if isfield(RUN,'run') && isfield(RUN.run, phase)
        t0_ts = RUN.run.(phase).t0_ts;
        t1_ts = RUN.run.(phase).t1_ts;
    end

    Units = [];
    if exist(uFile,'file')==2
        S = load(uFile,'Units');
        Units = S.Units;
    end

    % identify whether we have Baphy/epochs for this phase
    syncFile = fullfile(datapath, ['Master_sync_' phase '.mat']);
    hasHF = exist(syncFile,'file')==2;

    B = []; E = [];
    if hasHF
        Ssyn = load(syncFile);
        B = pick_struct_field(Ssyn, {'Baphy','BC','BP'});
        E = pick_struct_field(Ssyn, {'Epochs','EC','EP'});
    end

    % choose representative unit
    repIdx = NaN;
    zetaP = nan(numel(Units),1);

    if hasHF && ~isempty(Units) && isstruct(E) && isfield(E,'stim_start_ts') && isfield(E,'stim_stop_ts')
        stimStart_ts = E.stim_start_ts(:);
        stimStop_ts  = E.stim_stop_ts(:);
        good = isfinite(stimStart_ts) & isfinite(stimStop_ts) & (stimStop_ts > stimStart_ts);

        t0 = double(stimStart_ts(good))/TsRate;
        t1 = double(stimStop_ts(good))/TsRate;
        matEventTimes = [t0 t1];

        for u = 1:numel(Units)
            try
                zetaP(u) = zetatest(Units(u).spike_s(:), matEventTimes, 3.5, 500);
            catch
                zetaP(u) = NaN;
            end
        end

        [~,repIdx] = min(zetaP);
        if ~isfinite(repIdx), repIdx = 1; end
    elseif ~isempty(Units)
        % fallback: pick max FR within epoch (if bounds known), else first
        repIdx = 1;
    end

    QC.(phase).nUnits = numel(Units);
    QC.(phase).repIdx = repIdx;
    QC.(phase).zetaP = zetaP;

    % ---------------- Row 1/2: Trial-aligned raster + PSTH (rep unit) ----------------
    axR = subplot(nRows, nCols, (1-1)*nCols + c); 
    axP = subplot(nRows, nCols, (2-1)*nCols + c);

    if hasHF && ~isempty(Units) && isfinite(repIdx) && isfield(E,'stim_start_ts') && isfield(E,'stim_stop_ts')
        spk_ts = Units(repIdx).spike_ts(:);
        trig_ts = E.stim_start_ts(:);

        win_s = [-2 5];
        tStart_ts = round(win_s(1)*TsRate);
        tEnd_ts   = round(win_s(2)*TsRate);
        bin_ts    = round(0.02*TsRate);

        FMA_plot_rasterpeth_axes(axR, axP, spk_ts, trig_ts, tStart_ts, tEnd_ts, bin_ts);
        title(axR, [phase ' | rep ' num2str(repIdx)], 'Interpreter','none');
        ylabel(axR,'trials');
        ylabel(axP,'FR');
        xlabel(axP,'ms');
    else
        axes(axR); axis off; text(0.1,0.5,'no trials','FontSize',10);
        axes(axP); axis off;
        title(axR, phase, 'Interpreter','none');
    end

    % ---------------- Row 3: Population raster across time ----------------
    axPopR = subplot(nRows, nCols, (3-1)*nCols + c); 
    if ~isempty(Units) && isfinite(t0_ts)
        A = FMA_units_to_ts_rel(Units, t0_ts, t1_ts, 120); % first 120s for readability
        if ~isempty(A)
            RasterPlot(tsdArray(A), 'AxHandle', axPopR, 'LineWidth', 1);
            xlabel(axPopR,'time (ms, from epoch start)');
            ylabel(axPopR,'units');
        else
            axes(axPopR); axis off; text(0.1,0.5,'no spikes','FontSize',10);
        end
    else
        axes(axPopR); axis off; text(0.1,0.5,'no units / no bounds','FontSize',10);
    end

    % ---------------- Row 4: Binned raster heatmap (MakeQfromS) ----------------
    axH = subplot(nRows, nCols, (4-1)*nCols + c); 
    if ~isempty(Units) && isfinite(t0_ts)
        if contains(lower(phase),'sleep'), bin_s = 10; else, bin_s = 1; end
        [D, t_min] = FMA_binned_heatmap(Units, t0_ts, t1_ts, bin_s);
        if ~isempty(D)
            imagesc(axH, t_min, 1:size(D,2), D'); axis(axH,'xy'); colormap('redblue'), caxis([-3 3])
            xlabel(axH,'time (min)');
            ylabel(axH,'units');
        else
            axes(axH); axis off; text(0.1,0.5,'heatmap empty','FontSize',10);
        end
    else
        axes(axH); axis off; text(0.1,0.5,'no units / no bounds','FontSize',10);
    end

    % ---------------- Row 5: "STRF-like": sound_name x time PSTH map ----------------
%     axSTRF = subplot(nRows, nCols, (5-1)*nCols + c); 
%     if hasHF && ~isempty(Units) && isfinite(repIdx) && isstruct(B) && isfield(B,'trial') && isfield(B.trial,'sound_name') && isfield(E,'stim_start_ts')
%         try
%             FMA_soundname_psth_map(axSTRF, Units(repIdx).spike_ts(:), E.stim_start_ts(:), string(B.trial.sound_name(:)));
%             title(axSTRF,'sound x time map');
%             xlabel(axSTRF,'time (ms)');
%             ylabel(axSTRF,'sound');
%         catch
%             axes(axSTRF); axis off; text(0.1,0.5,'sound map failed','FontSize',10);
%         end
%     else
%         axes(axSTRF); axis off; text(0.1,0.5,'no sound_name','FontSize',10);
%     end

    % ---------------- Row 6: ZETA summary ----------------
    axZ = subplot(nRows, nCols, (6-1)*nCols + c); 
    if hasHF && ~isempty(zetaP)
        zp = zetaP(isfinite(zetaP));
        if isempty(zp)
            axes(axZ); axis off; text(0.1,0.5,'no zeta','FontSize',10);
        else
            histogram(axZ, zp, 30);
            hold(axZ,'on');
            yl = ylim(axZ);
            plot(axZ, [0.05 0.05], yl, '--k');
            hold(axZ,'off');
            frac = mean(zp<=0.05);
            title(axZ, sprintf('zeta frac<=0.05: %.2f', frac));
            xlabel(axZ,'p'); ylabel(axZ,'count');
        end
    else
        axes(axZ); axis off; text(0.1,0.5,'no zeta','FontSize',10);
    end

end

sgtitle(fh, ['FMA sanity | ' datapath], 'Interpreter','none');

pngOut = fullfile(outDir, 'FMA_epoch_sanity.png');
saveas(fh, pngOut);
save(fullfile(outDir, 'FMA_epoch_sanity_QC.mat'), 'QC', '-v7.3');
close(fh);

end

% ---------------- subfunctions ----------------

function out = pick_struct_field(S, names)
out = [];
for i = 1:numel(names)
    nm = names{i};
    if isfield(S,nm)
        out = S.(nm);
        return
    end
end
end

function FMA_plot_rasterpeth_axes(axRaster, axPsth, spk_ts, trig_ts, tStart_ts, tEnd_ts, bin_ts)
S = ts(spk_ts);
center = ts(trig_ts);

is = intervalSet(Range(center)+tStart_ts, Range(center)+tEnd_ts);
sweeps = intervalSplit(S, is, 'OffsetStart', tStart_ts);

ss = oneSeries(sweeps);
sq = intervalRate(ss, regular_interval(tStart_ts, tEnd_ts, bin_ts));

axes(axRaster);
cla(axRaster);
RasterPlot(sweeps, 'AxHandle', axRaster, 'TStart', tStart_ts, 'TEnd', tEnd_ts, 'LineWidth', 1);
set(axRaster,'Box','on');

axes(axPsth);
cla(axPsth);
dArea = Data(sq)/length(sweeps);
area(Range(sq,'ms'), dArea, 'FaceColor','k');
xlim(axPsth, [tStart_ts tEnd_ts]/10);
if max(dArea)>0, ylim(axPsth, [0 max(dArea)*1.2]); end
set(axPsth,'Box','on');
end

function A = FMA_units_to_ts_rel(Units, t0_ts, t1_ts, max_s)
% returns cell array of ts with times relative to t0_ts, restricted to [0..max_s]
TsRate = 1e4;
if nargin < 4, max_s = 120; end
tMax = round(max_s*TsRate);

A = {};
for u = 1:numel(Units)
    x = Units(u).spike_ts(:) - t0_ts;
    x = x(x>=0 & x<=min(tMax, (t1_ts-t0_ts)));
    if isempty(x), continue; end
    A{end+1} = ts(x); 
end
end

function [D, t_min] = FMA_binned_heatmap(Units, t0_ts, t1_ts, bin_s)
% MakeQfromS, zscore per unit
TsRate = 1e4;
bin_ts = round(bin_s*TsRate);
if bin_ts < 1, bin_ts = 1; end

A = {};
for u = 1:numel(Units)
    x = Units(u).spike_ts(:) - t0_ts;
    x = x(x>=0 & x<=(t1_ts-t0_ts));
    A{end+1} = ts(x); 
end
if isempty(A)
    D = []; t_min = []; return
end

B = tsdArray(A);
Q = MakeQfromS(B, bin_ts); % counts
D0 = full(Data(Q));        % [T x N]

% zscore per unit (column-wise), NaN-safe
mu = nanmean(D0,1);
sd = nanstd(D0,[],1);
sd(sd==0) = NaN;
D = (D0 - repmat(mu,size(D0,1),1)) ./ repmat(sd,size(D0,1),1);

t_min = (Range(Q)/TsRate)/60; % minutes
end

function FMA_soundname_psth_map(ax, spk_ts, trig_ts, sound_name)
% Build PSTH per unique sound_name (mean FR in bins), plot as image (sound x time)
TsRate = 1e4;
win_s = [-0.2 1.0];
bin_s = 0.02;

edges = (win_s(1):bin_s:win_s(2));
t_ms  = (edges(1:end-1)+bin_s/2)*1000;

uSounds = unique(sound_name);
uSounds(uSounds=="") = [];
if isempty(uSounds)
    axes(ax); axis off; text(0.1,0.5,'no sound names'); return
end

spk_s = double(spk_ts)/TsRate;
trig_s = double(trig_ts)/TsRate;

M = nan(numel(uSounds), numel(t_ms));

for iS = 1:numel(uSounds)
    idx = (sound_name == uSounds(iS));
    if ~any(idx), continue; end
    st = trig_s(idx);

    cnt = zeros(1,numel(edges)-1);
    for k = 1:numel(st)
        x = spk_s - st(k);
        cnt = cnt + histcounts(x, edges);
    end
    M(iS,:) = cnt / (numel(st)*bin_s); % Hz
end

imagesc(ax, t_ms, 1:numel(uSounds), M);
axis(ax,'xy');
end