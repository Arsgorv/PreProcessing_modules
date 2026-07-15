
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEAD RESTRAINT, 6 markers, small timescale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

cd('/media/nas8/OB_ferret_AG_BM/Shropshire/head-fixed/20241210_TORCs') % HR
% cd('/media/nas8/OB_ferret_AG_BM/Shropshire/freely-moving/20241205_TORCs') % FM

%% load variables
load('behavResources.mat', 'MovAcctsd')
load('HeartBeatInfo.mat', 'EKG')
load('SleepScoring_OBGamma.mat', 'SmoothGamma','Mean_FR','EMG_tsd', 'Clean_pupil_size')

Var{1} = MovAcctsd;
Var{2} = EKG.HBRate;
Var{3} = SmoothGamma;
Var{4} = Clean_pupil_size;
Var{5} = Mean_FR;
Var{6} = EMG_tsd;

Params = {'Motion','Heart rate','OB gamma power','Pupil size','Mean FR','EMG'};


%% clean data
D = Data(Var{1}); D([1 2]) = D(3); Var{1} = tsd(Range(Var{1}) , D); % correction first data accelero


%% parameters
set(0,'DefaultFigureWindowStyle','docked')

smootime_large = 1000; % in seconds
win_zsc = 1000; % in seconds
fs_plot = 10; % sampling for all
maxLagSec = 6; % for cross-corr in seconds


%% corr
for i=3:4%1:length(Params)
    
    Var{i} = Restrict(Var{i} , Var{5}); % every one binned at 10Hz;
    D = Data(Var{i});
    if or(or(i==1 , i==3) , i==6)
        D = log10(D);
    elseif i==4
        D = D/1e3;
    end
    D_small = D;
    
    Var_smooth_small{i} = tsd(Range(Var{i}) , zscore_sliding(D_small , win_zsc*fs_plot)); % detrending
    plot(Range(Var_smooth_small{i} , 's') , Data(Var_smooth_small{i})), hold on
end


% gathering data in X
clear X
nVar = numel(Params);
for i = 1:nVar
%     X(:,i) = zscore_nan(Data(Var_smooth_large{i}));
    X(:,i) = Data(Var_smooth_small{i});
end
X([1:20e3 end-5e3:end],:) = [];


% Correlation matrix
R = corr(X , 'Rows','pairwise');
R2 = R-R(3,4);
figure('Color','w');
subplot(131)
imagesc(R2); axis square; c=colorbar; ylabel(c,'R norm. to OB/pupil'); caxis([-.5 .5]);
set(gca,'XTick',1:nVar,'XTickLabel',Params,'YTick',1:nVar,'YTickLabel',Params);
xtickangle(45)
title('Pairwise correlation')
colormap redblue


% Cross-correlation vs reference
ref_name  = 'pupil';
idx_ref = find(contains(lower(Params), lower(ref_name)),1);
if isempty(idx_ref), idx_ref = 1; end
maxLag = round(maxLagSec * fs_plot);

X_slid_Zsc = zscore_sliding(X , win_zsc*10);

peakLagSec = zeros(1,nVar);
peakR = zeros(1,nVar);
for i = 1:nVar
    xi = X_slid_Zsc(:,i);
    xr = X_slid_Zsc(:,idx_ref);
    mask = isfinite(xi) & isfinite(xr);
    xi(~mask) = 0;
    xr(~mask) = 0;
    [xc,lags] = xcorr(xi, xr, maxLag, 'coeff');
    [~,ix] = max(abs(xc));
    peakLagSec(i) = lags(ix) / fs_plot;
    peakR(i) = xc(ix);
end

subplot(132)
bar(peakLagSec([1 2 3 5 6]) , 'FaceColor' , 'k')
set(gca,'XTick',1:nVar-1,'XTickLabel',Params([1 2 3 5 6])); xtickangle(45)
ylabel('Lag (s)')
title(['Peak lag vs ',Params{idx_ref}])
axis square, box off

subplot(133)
bar(peakR([1 2 3 5 6]) , 'FaceColor' , 'k')
set(gca,'XTick',1:nVar-1,'XTickLabel',Params([1 2 3 5 6])); xtickangle(45)
ylabel('Peak |R|')
axis square, box off



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEAD RESTRAINT, 6 markers, big timescale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

cd('/media/nas8/OB_ferret_AG_BM/Shropshire/head-fixed/20241210_TORCs') % HR
% cd('/media/nas8/OB_ferret_AG_BM/Shropshire/freely-moving/20241205_TORCs') % FM

%% load variables
load('behavResources.mat', 'MovAcctsd')
load('HeartBeatInfo.mat', 'EKG')
load('SleepScoring_OBGamma.mat', 'SmoothGamma','Mean_FR','EMG_tsd', 'Clean_pupil_size')

Var{1} = MovAcctsd;
Var{2} = EKG.HBRate;
Var{3} = SmoothGamma;
Var{4} = Clean_pupil_size;
Var{5} = Mean_FR;
Var{6} = EMG_tsd;

Params = {'Motion','Heart rate','OB gamma power','Pupil size','Mean FR','EMG'};


%% clean data
D = Data(Var{1}); D([1 2]) = D(3); Var{1} = tsd(Range(Var{1}) , D); % correction first data accelero


%% parameters
set(0,'DefaultFigureWindowStyle','docked')

smootime_small = 1; % in seconds
smootime_large = 1000; % in seconds
win_zsc = 1000; % in seconds
fs_plot = 10; % sampling for all
maxLagSec = 180; % for cross-corr in seconds


%% corr
for i=1:length(Params)
    
    Var{i} = Restrict(Var{i} , Var{5}); % every one binned at 10Hz;
    D = Data(Var{i});
    if or(or(i==1 , i==3) , i==6)
        D = log10(D);
    end
    D_large = movmean(D , smootime_large*fs_plot , 'omitnan');
%     D_small = movmean(D , smootime_small*10 , 'omitnan');
    D_small = D;
    
    Var_smooth_large{i} = tsd(Range(Var{i}) , D_large);
    Var_smooth_small{i} = tsd(Range(Var{i}) , zscore_sliding(D_small , win_zsc*fs_plot));
end


% gathering data in X
clear X
nVar = numel(Params);
for i = 1:nVar
    X(:,i) = zscore_nan(Data(Var_smooth_large{i}));
%     X(:,i) = Data(Var_smooth_small{i});
end
X([1:20e3 end-5e3:end],:) = [];


%%
% --- Parameters ---
ref_name  = 'pupil';
idx_ref = find(contains(lower(Params), lower(ref_name)),1);
if isempty(idx_ref), idx_ref = 1; end
maxLag = round(maxLagSec * fs_plot);

% Sliding z-score
X_slid_Zsc = zscore_sliding(X , win_zsc*10);

nVar = size(X_slid_Zsc,2);
peakLagSec = nan(nVar,nVar);
peakR = nan(nVar,nVar);

% --- Loop over all pairs ---
for i = 1:nVar
    for j = 1:nVar
        xi = X_slid_Zsc(:,i);
        xj = X_slid_Zsc(:,j);
        mask = isfinite(xi) & isfinite(xj);
        xi(~mask) = 0;
        xj(~mask) = 0;

        [xc,lags] = xcorr(xi, xj, maxLag, 'coeff');
        [~,ix] = max(abs(xc));
        peakLagSec(i,j) = lags(ix) / fs_plot;  % lag in seconds
        peakR(i,j)      = xc(ix);              % correlation value
    end
end

% --- Plot ---
figure('Color','w');

% Subplot 1: Peak lag matrix
subplot(1,2,1)
[Data_corr3 , ~ , ~ , v1] = OrderMatrix_BM(peakLagSec  , Params , [] , 0);
imagesc(Data_corr3); axis square
c=colorbar; ylabel(c,'Lag (s)');
set(gca,'XTick',1:nVar,'XTickLabel',Params(v1), ...
        'YTick',1:nVar,'YTickLabel',Params(v1));
xtickangle(45)
title('Peak lag (s)')
colormap redblue

% Subplot 2: Peak |R| matrix
subplot(1,2,2)
peakR2 = peakR-peakR(3,4);
imagesc(peakR2(v1,v1)); axis square
c=colorbar; ylabel(c,'Peak |R|');
set(gca,'XTick',1:nVar,'XTickLabel',Params(v1), ...
        'YTick',1:nVar,'YTickLabel',Params(v1));
xtickangle(45), caxis([-.5 .5])
title('Peak correlation')
colormap redblue





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FREELY MOVING, 5 markers, -HR/pupil +sleep states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;
set(0,'DefaultFigureWindowStyle','docked')

% cd('/media/nas8/OB_ferret_AG_BM/Shropshire/head-fixed/20241210_TORCs') % HR
cd('/media/nas8/OB_ferret_AG_BM/Shropshire/freely-moving/20241206_TORCs') % FM

%% --- Load variables ---
load('behavResources.mat', 'MovAcctsd')
load('SleepScoring_OBGamma.mat', ...
    'SmoothGamma_wide','Mean_FR','EMG_tsd','SmoothTheta','SmoothDelta_OB','Wake','Sleep')

% TEMP redefine epochs if needed
% Wake  = or(intervalSet(0,1000*1e4), intervalSet(14000*1e4,15200*1e4));
Sleep = and(Sleep , intervalSet(5000*1e4, 14000*1e4));

Var = {MovAcctsd, SmoothGamma_wide, Mean_FR, EMG_tsd, SmoothTheta , SmoothDelta_OB};
Params = {'Motion','OB gamma power','Mean FR','EMG','HPC theta power','OB delta power'};

%% --- Parameters ---
win_zsc_list  = [500, 10000];     % z-score windows (samples)
win_sm_list   = [50, 300];        % smoothing windows (samples)
xlims_list    = {[-300 300], [-2400 2400]};
maxLag_list   = [3000, 24000];    % max lags (samples)
fs_plot       = 10;               % sampling rate (Hz)
nVar          = numel(Var);

conds = {'All','Wake','Sleep'};
cond_epochs = {[], Wake, Sleep};  % empty for "All"

%% --- Reference variable ---
ref_idx = find(strcmpi(Params, 'Mean FR'), 1);
if isempty(ref_idx), error('Mean FR not found in Params'); end

%% --- Storage for peaks ---
peakLagSec = cell(numel(win_zsc_list), numel(conds));
peakR      = cell(numel(win_zsc_list), numel(conds));

%% --- Loop over z-score windows ---
for w = 2%1:numel(win_zsc_list)
    win_zsc = win_zsc_list(w);
    win_sm  = win_sm_list(w);
    maxLag  = maxLag_list(w);
    xlims   = xlims_list{w};

    %% --- Loop over conditions ---
    for c = 3%1:numel(conds)

        % Restrict variables to current state
        if c==1
            Ep = intervalSet(0 , 15.5e7);
        elseif c==2
            Ep = Wake;
        elseif c==3
            Ep = Sleep;
        end
        for i = 1:nVar
            Var_state{i} = Restrict(Var{i}, Ep);
        end
        for i = 1:nVar
            Var_state{i} = Restrict(Var_state{i}, Var_state{3});
        end
        
        % Smooth after restriction
        X_smooth = nan(length(Data(Var_state{i})), nVar);
        for i = 1:nVar
            if or(or(i==1 , i==2) , i==4)
                r = movmean(log10(Data(Var_state{i})), win_sm, 'omitnan')';
            else
                r = movmean(Data(Var_state{i}), win_sm, 'omitnan')';
            end
            r(r == -Inf) = NaN;
            X_smooth(:,i) = r;
        end

        % Z-score with chosen window
        Xc = nan(size(X_smooth));
        for i = 1:nVar
            Xc(:,i) = zscore_sliding(X_smooth(:,i), win_zsc);
        end

        %% --- Cross-correlation ---
        peakLagSec{w,c} = nan(nVar);
        peakR{w,c}      = nan(nVar);
%         figure('Name',sprintf('Cross-corr (%s, win=%d)',conds{c},win_zsc),'Color','w')
        for i = 1:nVar
            for j = 2%i+1:nVar
                xi = Xc(:,i); xj = Xc(:,j);
                mask = isfinite(xi) & isfinite(xj);
                xi(~mask) = 0; xj(~mask) = 0;
%                 [xc,lags] = xcorr(xi, xj, maxLag, 'coeff');
                [xc(i,:),lags] = xcorr(xi, xj, maxLag, 'coeff');
%                 [~,ix] = max(abs(xc));
%                 peakLagSec{w,c}(i,j) = lags(ix)/fs_plot;
%                 peakR{w,c}(i,j)      = xc(ix);
%                 subplot(nVar,nVar,(i-1)*nVar + j)
%                 plot(lags/fs_plot, xc, 'k')
%                 xlabel('Lag (s)'); ylabel('R')
%                 xlim(xlims); grid on
%                 title(sprintf('%s vs %s',Params{i},Params{j}))
            end
        end

        %% --- Autocorrelation ---
%         figure('Name',sprintf('Autocorr (%s, win=%d)',conds{c},win_zsc),'Color','w')
%         for i = 1:nVar
%             xi = Xc(:,i); xi(~isfinite(xi)) = 0;
%             [ac,lags] = xcorr(xi, xi, maxLag, 'coeff');
%             subplot(1,nVar,i)
%             plot(lags/fs_plot, ac, 'k')
%             xlim(xlims); grid on
%             title(Params{i})
%             xlabel('Lag (s)'); ylabel('R')
%         end
    end
end

figure, hold on
for i=[1 3 4 5 6]
    plot(lags/fs_plot, xc(i,:))
    xlabel('Lag (s)'); ylabel('R')
    xlim(xlims); grid on
end
legend(Params([1 3 4 5 6]))
ylim([-1 1])




%% nice display for 10/12 HR

figure, hold on
% subplot(121)
for i=1:length(Params)
    
    Var{i} = Restrict(Var{i} , Var{5});
    Var{i} = Restrict(Var{i} , intervalSet(3200e4 , 4800e4));
    D = Data(Var{i});
    
    if or(or(i==1 , i==3) , i==6)
        D = log10(D);
    end
    
    D = [movmean(D(1:7030) ,  100*10 , 'omitnan') ;...
        movmean(D(7031:8030) , 20*10 , 'omitnan');...
        movmean(D(8031:end) , 100*10 , 'omitnan')];
    D = zscore_nan(D);
    
    if or(i==3 , i==4)
        plot(linspace(0,1600,length(D))/60 , D+1)
    else
        plot(linspace(0,1600,length(D))/60 , D)
    end
end
legend(Params)
xlabel('Time (min)'), ylabel('Norm. values (a.u.'), xlim([0 25]), ylim([-5 7])
box off






%% nice display for freely moving
clear all

% cd('/media/nas8/OB_ferret_AG_BM/Shropshire/head-fixed/20241210_TORCs') % HR
cd('/media/nas8/OB_ferret_AG_BM/Shropshire/freely-moving/20241205_TORCs') % FM

%% load variables
load('behavResources.mat', 'MovAcctsd')
load('SleepScoring_OBGamma.mat', 'SmoothGamma','Mean_FR','EMG_tsd','SmoothTheta' , 'Wake', 'Sleep')

Var{1} = MovAcctsd;
Var{2} = SmoothGamma;
Var{3} = Mean_FR;
Var{4} = EMG_tsd;
Var{5} = SmoothTheta;

Params = {'Motion','OB gamma power','Mean FR','EMG','HPC theta power'};

%% --- Parameters ---
win_sm  = 200;      % smoothing window (samples)
win_zsc = 10000;    % z-score window (samples)
maxLagSec = 600;    % max lag in seconds for cross-correlation/autocorrelation
fs_plot = 10;       % assumed sampling rate (Hz)
nVar = numel(Var);


%% --- Find reference variable index ---
ref_idx = find(strcmpi(Params, 'Mean FR'), 1);
if isempty(ref_idx)
    error('Mean FR not found in Params');
end

% --- Restrict all variables to reference timestamps ---
for i = 1:nVar
    Var{i} = Restrict(Var{i}, Var{ref_idx});
end
t_ref = Range(Var{ref_idx}, 's');

% --- Smooth + sliding z-score for all variables ---
X_all = nan(numel(t_ref), nVar);
for i = 1:nVar
    r = movmean(log10(Data(Var{i})), win_sm, 'omitnan')';
    r(r == -Inf) = NaN;
    X_all(:,i) = zscore_sliding(r, win_zsc);
end

% --- Prepare conditions ---
conds = {'All','Wake','Sleep'};
X_cond = cell(1, numel(conds));
t_cond = cell(1, numel(conds));

% All
X_cond{1} = X_all;
t_cond{1} = t_ref;

% Wake
Wake_tsd = Restrict(Var{ref_idx}, Wake);
t_cond{2} = Range(Wake_tsd, 's');
[~, idx_wake] = ismember(t_cond{2}, t_ref);
idx_wake(idx_wake == 0) = [];
X_cond{2} = X_all(idx_wake, :);

% Sleep
Sleep_tsd = Restrict(Var{ref_idx}, Sleep);
t_cond{3} = Range(Sleep_tsd, 's');
[~, idx_sleep] = ismember(t_cond{3}, t_ref);
idx_sleep(idx_sleep == 0) = [];
X_cond{3} = X_all(idx_sleep, :);

% =========================
% --- Plot time series ---
% =========================
figure('Color','w','Name','Z-scored variables over time');
for c = 1:numel(conds)
    subplot(3,1,c); hold on
    for i = 1:nVar
        plot(t_cond{c}, X_cond{c}(:,i), 'LineWidth', 1);
    end
    title(conds{c});
    ylabel('z-score (a.u.)');
    if c == 3
        xlabel('Time (s)');
    end
    legend(Params, 'Interpreter','none');
    grid on
end
