function PlotSessionMetadata_OB(sessions, save_dir)
% PlotSessionMetadata_OB  Visualize per-session metadata for OB ferret recordings.
%
% GOAL
%   For each animal, produce a single figure of control panels for the
%   paper, focused on how recordings are distributed across the day
%   (circadian consistency) and on basic animal-state covariates
%   (weight, age, stabulation temperature). Also a cross-animal
%   comparison figure.
%
% USAGE
%   PlotSessionMetadata_OB                              % default sessions, save to pwd
%   PlotSessionMetadata_OB(sessions)                    % sessions = cell of paths
%   PlotSessionMetadata_OB(Dir.path, 'C:\figs')         % save to folder
%   PlotSessionMetadata_OB(Dir, 'C:\figs')              % accepts PathForExperimentsOB output
%
% INPUT
%   sessions  - cell array of session folder paths, OR a Dir struct with
%               a .path cell. If empty, builds the default list:
%               PathForExperimentsOB({'Brynza','Labneh','Shropshire'},
%               {'freely-moving','head-fixed'}).
%   save_dir  - folder where figures are saved (.fig + .png).
%               Default: pwd.
%
% OUTPUT
%   One figure per animal + one cross-animal figure. Saved as
%   SessionMetadata_<animal>.{fig,png} and SessionMetadata_comparison.{fig,png}.
%
% Arsenii / OB project.

%% -------------------- input handling --------------------
if nargin < 1 || isempty(sessions)
    D1 = PathForExperimentsOB({'Brynza','Labneh','Shropshire'}, 'freely-moving');
    D2 = PathForExperimentsOB({'Brynza','Labneh','Shropshire'}, 'head-fixed');
    sessions = unique([D1.path, D2.path], 'stable');
end
if isstruct(sessions) && isfield(sessions, 'path')
    sessions = sessions.path;
end
if ~iscell(sessions), sessions = {sessions}; end

if nargin < 2 || isempty(save_dir), save_dir = pwd; end
if exist(save_dir, 'dir') ~= 7, mkdir(save_dir); end

%% -------------------- load metadata --------------------
M = collect_metadata(sessions);
if isempty(M)
    error('No session_metadata.mat files found in the provided sessions.');
end
fprintf('Loaded metadata from %d sessions.\n', numel(M));

%% -------------------- per-animal figures --------------------
animals = unique({M.animal});
for a = 1:numel(animals)
    aname = animals{a};
    Ma = M(strcmp({M.animal}, aname));
    if isempty(Ma), continue; end
    plot_one_animal(Ma, aname, save_dir);
end

%% -------------------- cross-animal comparison --------------------
plot_comparison(M, animals, save_dir);

fprintf('Figures saved to %s\n', save_dir);
end


%% ===================================================================
%% Loading
%% ===================================================================

function M = collect_metadata(sessions)
% Project each session_metadata.mat onto a fixed schema so struct
% concatenation is robust to extra/missing fields between files.
fields  = {'animal','animal_age_days','animal_weight_g','mean_temperature_c', ...
           'recording_date','start_time','end_time','duration_sec','recording_path'};
strFlds = {'animal','recording_date','recording_path'};
dtFlds  = {'start_time','end_time'};

M = struct([]);
for i = 1:numel(sessions)
    f = fullfile(sessions{i}, 'session_metadata.mat');
    if exist(f, 'file') ~= 2
        warning('Skipping (no session_metadata.mat): %s', sessions{i});
        continue;
    end
    S = load(f, 'session_metadata');
    if ~isfield(S, 'session_metadata'), continue; end
    sm = S.session_metadata;
    out = struct();
    for k = 1:numel(fields)
        fn = fields{k};
        if isfield(sm, fn)
            out.(fn) = sm.(fn);
        elseif any(strcmp(fn, dtFlds))
            out.(fn) = NaT;
        elseif any(strcmp(fn, strFlds))
            out.(fn) = '';
        else
            out.(fn) = NaN;
        end
    end
    if isempty(M), M = out; else, M(end+1) = out; end %#ok<AGROW>
end
end


%% ===================================================================
%% Per-animal figure
%% ===================================================================

function plot_one_animal(M, aname, save_dir)
% 2x3 panels: weight+age timeline | circadian polar | date x hour heatmap |
%             duration hist | temperature timeline | correlation matrix

% --- extract numeric arrays, dropping NaT/NaN where needed ---
nS = numel(M);
dnum   = nan(nS, 1);
sHour  = nan(nS, 1);
eHour  = nan(nS, 1);
durH   = nan(nS, 1);
wt     = nan(nS, 1);
age    = nan(nS, 1);
tempC  = nan(nS, 1);

for i = 1:nS
    if isa(M(i).start_time,'datetime') && ~isnat(M(i).start_time)
        dt = M(i).start_time;
        dnum(i)  = datenum(dt);
        sHour(i) = hour(dt) + minute(dt)/60 + second(dt)/3600;
    end
    if isa(M(i).end_time,'datetime') && ~isnat(M(i).end_time)
        dt = M(i).end_time;
        eHour(i) = hour(dt) + minute(dt)/60 + second(dt)/3600;
    end
    if ~isnan(M(i).duration_sec), durH(i) = M(i).duration_sec/3600; end
    if ~isnan(M(i).animal_weight_g),    wt(i)    = M(i).animal_weight_g; end
    if ~isnan(M(i).animal_age_days),    age(i)   = M(i).animal_age_days; end
    if ~isnan(M(i).mean_temperature_c), tempC(i) = M(i).mean_temperature_c; end
end

% sort by recording date
[~, ord] = sort(dnum);
dnum  = dnum(ord);  sHour = sHour(ord);  eHour = eHour(ord);
durH  = durH(ord);  wt    = wt(ord);     age   = age(ord);
tempC = tempC(ord);

% --- figure ---
hF = figure('Name', sprintf('SessionMetadata - %s', aname), ...
            'Color', 'w', 'Position', [60 60 1500 850]);

% ---- (1) weight + age timeline ----
ax1 = subplot(2, 3, 1);
yyaxis left
plot(dnum, wt, 'o-', 'LineWidth', 1.2, 'MarkerFaceColor', 'auto'); hold on;
ylabel('Body weight (g)');
yyaxis right
plot(dnum, age, 's-', 'LineWidth', 1.2);
ylabel('Age at recording (days)');
xlabel('Date');
title('Weight & age over time');
datetick(ax1, 'x', 'yyyy-mm-dd', 'keeplimits');
xtickangle(ax1, 30);
grid on;

% ---- (2) circadian polar histogram of start times ----
ax2 = subplot(2, 3, 2, polaraxes);
ok = ~isnan(sHour);
if any(ok)
    theta = 2*pi * sHour(ok) / 24;     % 0-24h -> 0-2pi
    polarhistogram(theta, 24, 'FaceColor', [0.25 0.45 0.75], ...
                   'EdgeColor', 'k', 'FaceAlpha', 0.8);
end
ax2.ThetaZeroLocation = 'top';
ax2.ThetaDir = 'clockwise';
ax2.ThetaTick = 0:45:315;
ax2.ThetaTickLabel = {'0h','3h','6h','9h','12h','15h','18h','21h'};
title(sprintf('Start time circadian (n=%d)', sum(ok)));

% ---- (3) date x hour heatmap of recording occupation ----
ax3 = subplot(2, 3, 3);
[H, dayLabels] = build_day_hour_matrix(dnum, sHour, eHour);
imagesc(0:23, 1:size(H,1), H);
colormap(ax3, custom_blue());
xlabel('Hour of day');
ylabel('Recording date');
set(ax3, 'YTick', 1:size(H,1), 'YTickLabel', dayLabels, 'FontSize', 8);
xlim([-0.5 23.5]);
title('Recording occupation (date \times hour)');
hb = colorbar; ylabel(hb, 'fraction of hour recorded');

% ---- (4) duration distribution ----
ax4 = subplot(2, 3, 4);
ok = ~isnan(durH);
if any(ok)
    histogram(durH(ok), max(5, round(sum(ok)/3)), 'FaceColor', [0.4 0.7 0.4]);
end
xlabel('Recording duration (h)');
ylabel('Count');
title(sprintf('Durations  median=%.2f h', nanmedian_local(durH)));
grid on;

% ---- (5) temperature timeline ----
ax5 = subplot(2, 3, 5);
ok = ~isnan(tempC);
if any(ok)
    plot(dnum(ok), tempC(ok), 'o-', 'Color', [0.8 0.4 0.2], 'LineWidth', 1.2, ...
         'MarkerFaceColor', [0.8 0.4 0.2]);
end
xlabel('Date');
ylabel('Mean stabulation T (\circC)');
title('Temperature across sessions');
datetick(ax5, 'x', 'yyyy-mm-dd', 'keeplimits');
xtickangle(ax5, 30);
grid on;

% ---- (6) correlation matrix ----
ax6 = subplot(2, 3, 6);
plot_correlation_matrix(ax6, [wt age sHour durH tempC], ...
    {'weight','age','start hr','dur (h)','temp'});

% ---- supertitle ----
suptitle_local(sprintf('%s    (%d sessions, %.1f h total)', ...
    aname, nS, nansum_local(durH)));

% ---- save ----
out_png = fullfile(save_dir, sprintf('SessionMetadata_%s.png', aname));
out_fig = fullfile(save_dir, sprintf('SessionMetadata_%s.fig', aname));
try, saveas(hF, out_png); catch, end
try, savefig(hF, out_fig); catch, end
fprintf('  %s: %d sessions -> %s\n', aname, nS, out_png);
end


%% ===================================================================
%% Cross-animal comparison
%% ===================================================================

function plot_comparison(M, animals, save_dir)
% 2x2 panels: n sessions | start-hour boxplot | weight boxplot | age boxplot

nA = numel(animals);
nS = zeros(nA, 1);
sH_all = []; wt_all = []; age_all = []; dur_all = []; grp = [];

for a = 1:nA
    Ma = M(strcmp({M.animal}, animals{a}));
    nS(a) = numel(Ma);
    for i = 1:numel(Ma)
        s = Ma(i);
        if isa(s.start_time,'datetime') && ~isnat(s.start_time)
            sH_all(end+1,1) = hour(s.start_time) + minute(s.start_time)/60; %#ok<AGROW>
        else
            sH_all(end+1,1) = NaN; %#ok<AGROW>
        end
        wt_all (end+1,1) = s.animal_weight_g;   %#ok<AGROW>
        age_all(end+1,1) = s.animal_age_days;   %#ok<AGROW>
        if ~isnan(s.duration_sec)
            dur_all(end+1,1) = s.duration_sec/3600;  %#ok<AGROW>
        else
            dur_all(end+1,1) = NaN; %#ok<AGROW>
        end
        grp(end+1,1) = a; %#ok<AGROW>
    end
end

hF = figure('Name', 'SessionMetadata - comparison', ...
            'Color', 'w', 'Position', [80 80 1300 700]);

subplot(2, 2, 1);
bar(1:nA, nS, 'FaceColor', [0.3 0.5 0.8]);
set(gca, 'XTick', 1:nA, 'XTickLabel', animals);
ylabel('# sessions'); title('Sessions per animal'); grid on;
for a = 1:nA
    text(a, nS(a), sprintf('%d', nS(a)), 'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom', 'FontWeight','bold');
end

subplot(2, 2, 2);
boxplot_safe(sH_all, grp, animals);
ylabel('Recording start (h of day)');
title('Start time of day');
ylim([0 24]); set(gca, 'YTick', 0:3:24); grid on;

subplot(2, 2, 3);
boxplot_safe(wt_all, grp, animals);
ylabel('Body weight (g)'); title('Weight at recording'); grid on;

subplot(2, 2, 4);
boxplot_safe(age_all, grp, animals);
ylabel('Age (days)'); title('Age at recording'); grid on;

suptitle_local('Cross-animal comparison');

out_png = fullfile(save_dir, 'SessionMetadata_comparison.png');
out_fig = fullfile(save_dir, 'SessionMetadata_comparison.fig');
try, saveas(hF, out_png); catch, end
try, savefig(hF, out_fig); catch, end
end


%% ===================================================================
%% Small helpers
%% ===================================================================

function [H, dayLabels] = build_day_hour_matrix(dnum, sHour, eHour)
% Build a (nDates x 24) matrix; each cell = fraction of that hour recorded.
ok = ~isnan(dnum);
days = floor(dnum(ok));
uDays = unique(days);
H = zeros(numel(uDays), 24);
for i = 1:numel(ok)
    if ~ok(i), continue; end
    d   = floor(dnum(i));
    row = find(uDays == d, 1);
    sH  = sHour(i);
    eH  = eHour(i);
    if isnan(sH) || isnan(eH), continue; end
    if eH < sH, eH = 24; end                  % treat as ending at midnight
    for h = 0:23
        a = max(sH, h);
        b = min(eH, h+1);
        if b > a, H(row, h+1) = b - a; end    % overlap of [sH,eH] with [h,h+1]
    end
end
dayLabels = cellstr(datestr(uDays, 'yyyy-mm-dd'));
end


function plot_correlation_matrix(ax, X, labels)
% Pairwise correlation of columns of X, rendered as a labeled heatmap.
axes(ax); %#ok<LAXES>
% guard: drop all-NaN columns
keep = ~all(isnan(X), 1);
X = X(:, keep);
labels = labels(keep);
if size(X, 2) < 2
    text(0.5, 0.5, 'not enough data', 'HorizontalAlignment','center');
    axis off; return;
end
try
    R = corr(X, 'rows', 'pairwise');
catch
    R = nan(size(X, 2));
    for i = 1:size(X,2)
        for j = 1:size(X,2)
            a = X(:,i); b = X(:,j);
            ok = ~isnan(a) & ~isnan(b);
            if sum(ok) >= 3
                c = corrcoef(a(ok), b(ok));
                R(i,j) = c(1,2);
            end
        end
    end
end
imagesc(R, [-1 1]);
colormap(ax, redblue_local());
hb = colorbar(ax); ylabel(hb, 'Pearson r');
n = size(R, 1);
set(ax, 'XTick', 1:n, 'YTick', 1:n, ...
        'XTickLabel', labels, 'YTickLabel', labels, ...
        'XTickLabelRotation', 30);
for i = 1:n
    for j = 1:n
        if isnan(R(i,j)), continue; end
        col = 'k'; if abs(R(i,j)) > 0.5, col = 'w'; end
        text(j, i, sprintf('%.2f', R(i,j)), 'HorizontalAlignment','center', ...
             'Color', col, 'FontSize', 9);
    end
end
title('Correlations');
axis square;
end


function boxplot_safe(y, grp, labels)
ok = ~isnan(y);
if ~any(ok), text(0.5, 0.5, 'no data', 'HorizontalAlignment','center'); axis off; return; end
boxplot(y(ok), grp(ok), 'Labels', labels);
end


function cm = redblue_local()
% Diverging red-white-blue colormap (64 rows).
n = 32;
r = [linspace(0.13, 1, n)'; linspace(1, 0.78, n)'];
g = [linspace(0.35, 1, n)'; linspace(1, 0.15, n)'];
b = [linspace(0.65, 1, n)'; linspace(1, 0.18, n)'];
cm = [r g b];
end


function cm = custom_blue()
% Sequential white-to-blue for the day x hour heatmap.
n = 64;
cm = [linspace(1, 0.18, n)', linspace(1, 0.35, n)', linspace(1, 0.65, n)'];
end


function suptitle_local(s)
% R2018b has no built-in suptitle in base MATLAB. Add a big invisible axes.
ax = axes('Position', [0 0.94 1 0.05], 'Visible', 'off');
text(0.5, 0.5, s, 'Parent', ax, 'HorizontalAlignment','center', ...
     'FontSize', 13, 'FontWeight','bold', 'Interpreter','none');
end


function m = nanmedian_local(x)
x = x(~isnan(x));
if isempty(x), m = NaN; else, m = median(x); end
end


function s = nansum_local(x)
x = x(~isnan(x));
if isempty(x), s = 0; else, s = sum(x); end
end