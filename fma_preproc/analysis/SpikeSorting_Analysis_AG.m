function SpikeSorting_Analysis_AG(directory)
% SpikeSorting_Analysis_AG
%
% A single script that:
%  1) Loads and preprocesses data from wave_clus or KlustaKwik, 
%  2) Performs spike rasters, PSTH analyses (raw & z-score), 
%  3) Generates population PSTHs, violin plots, heatmaps,
%  4) Explores correlations with LFP gamma power across vigilance states,
%  5) Runs control analyses (trial sub-sampling, random onsets, spike jitter),
%  6) Visualizes correlation matrices, box plots, etc.
%
% USAGE:
%   Master_SpikeAnalysis_FerretData('/path/to/sessionFolder');
%
% Options: 
% 20241204_TORCs, 20241205_TORCs, 20241206_TORCs, 20241209_TORCs, 20241212_TORCs, 20241213_TORCs, 20241214_TORCs, 20241221_TORCs, 20241223_TORCs
% 
% ToDo: 
%       -- Add Zeta-test to preselect spikes
%       -- Save raster plots to quickly access every channel
%       -- Add STRF calculation
% By: A. Goriachenkov, 01/2025
% -------------------------------------------------------------------------

%% SECTION 1.0: Basic Setup & Parameters
% if nargin < 1
%    error('Please provide a directory path to your ferret session data.');
% end
Dir = PathForExperimentsOB('Shropshire', 'freely-moving', 'all', 'TORCs');
Dir = PathForExperimentsOB('Shropshire', 'head-fixed', 'none', 'TORCs');

sessions = {...
            '20241205_TORCs', ...
            '20241206_TORCs', ...
            '20241209_TORCs', ...
            '20241210_TORCs', ...
            '20241211_TORCs', ...
            '20241212_TORCs', ...
            '20241213_TORCs', ...
            '20241214_TORCs', ...
            '20241221_TORCs', ...
            '20241223_TORCs'};

sessions = Dir.name';
exp_path  = '/media/nas8/OB_ferret_AG_BM/Shropshire/head-fixed/';
% '/media/nas8/OB_ferret_AG_BM/Shropshire/freely-moving/';
session = sessions{6}; %1

directory = [exp_path session];
cd(directory)
rasterDir = fullfile([directory '/wave_clus/'], 'raster_figures');
if ~exist(rasterDir, 'dir'), mkdir(rasterDir); end

%% IN PROGRESS:       SECTION 1.1: Apply ZETA test to select spiking channels
nSessions = length(sessions);
for s = 1:nSessions
    sessName = sessions{s};
    directory = [exp_path sessName];
    run_zeta_test_AG(directory, sessName)
end


%% SECTION 1.2: Load wave_clus Spikes (or KlustaKwik if you prefer)

[spikes, metadata] = load_wave_clus(directory);

% Optional: Remove bad data for a specific session
% if strcmp(session, '20241205_TORCs')
%     metadata([64:75], :) = [];
%     spikes(:, [64:75])   = [];
% end

%% SECTION 1.3: Load LFP Data & Extract Stimulus Onsets
lfpFile = fullfile(directory, 'LFPData', 'LFP111.mat');
if exist(lfpFile,'file')
    load(lfpFile,'LFP');
else
    error('LFP file not found: %s', lfpFile);
end

% Extract stimulus onsets:
StartSound = thresholdIntervals(LFP, 4000, 'Direction','Above');
Starttime  = Start(StartSound); 

% set(0, 'DefaultFigureWindowStyle','docked');  % or 'normal'

% Convert to convinient format for Yves (STRF)
% stim_start = Starttime/1e4;
% save([exp_path 'stim_starts/' 'stim_start_' session], 'stim_start')

%% FIGURE: Quick Raster Plots for Each Cluster

% Below is a loop that draws PSTHs for each cluster, to quickly check data:
Resp = [];  % store PSTH data if desired
for ii = 1:size(spikes,2)
    fh = figure('Visible', 'off');
    [fh,sq,sweeps] = RasterPETH(ts(spikes(:,ii)*1e1), ts(Starttime), ...
        -20000, +100000, 'BinSize', 240);
    i = metadata(ii,:);
    title(sprintf('Channel #%d, Cluster #%d', i(1), i(2)));
    if ~isempty(sq)
        Resp(ii,:) = Data(sq);
    end
    % Save figure for quick access
    fName = sprintf('Raster_ch%d_cl%d.png', i(1), i(2));
%     saveas(fh, fullfile(rasterDir, fName));
%     close(fh); 
end

%% Selecting â€œGoodâ€? Clusters (colIndex) Based on zeta test

clear chcl
chcl = get_chcl(directory, session);  % function that returns list of [channel cluster]

colIndex = [];
for s = 1:length(chcl)
    % Find row in metadata matching [channel cluster]
    row = find(metadata(:,1)==chcl{s}(1) & metadata(:,2)==chcl{s}(2));
    colIndex(s) = row;
end
% spikeTimes = spikes(:, colIndex);
% spikeTimes = spikeTimes(~isnan(spikeTimes)); % Remove NaN values

%% SECTION 1.5: Create a Cell Array A{n} for Neurons
% A{n} = ts(...) for each selected cluster
A = {};
n=1;
for i = 1:47;%colIndex
    A{n} = ts(spikes(:,i)*1e1);
    n = n+1;
end

%% FIGURE: quick raster of all chosen clusters
figure('Name','Raster of Chosen Clusters');
RasterPlot(A);

%% SECTION 1.6: Load Sleep & Gamma Data
% Load SleepScoring info, e.g. Epoch, Wake, Sleep, SWSEpoch, REMEpoch, etc.
sleepFile = fullfile(directory,'SleepScoring_OBGamma.mat');
if exist(sleepFile,'file')
    load(sleepFile,'Wake','Sleep','SWSEpoch','REMEpoch','Epoch_S1','Epoch_S2','SmoothGamma','Epoch');
else
    warning('No SleepScoring_OBGamma.mat found in %s', directory);
end

%% IN PROGRESSS:     SECTION 2: STRF calculation
STRF_calculation(directory, Starttime/1e4)


%% SECTION 3.1: Binning of spikes
clear B Q D

Binsize = 1*1e4;  % or 120*1e4
B = tsdArray(A);
Q = MakeQfromS(B,Binsize);
Q = tsd(Range(Q),nanzscore(full(Data(Q))));
D = Data(Q);

%% FIGURE: Binned raster heatmap
clear fh
fh = figure('Name','Binned Raster Heatmap', 'Units','normalized','Position',[0.1 0.3 0.7 0.4]);
subplot(3,1,1:2)
imagesc([1:size(D,1)]*(Binsize/1e4)/60, 1:size(D,2), D');
axis xy; caxis([-2 2]); colormap viridis
ylabel('clusters');
% c = colorbar; c.Label.String='firing rate (zscore)';
makepretty

subplot(3,1,3)
plot(linspace(0,size(D,1)*Binsize/1e4/60,size(D,1)), movmean(nanmean(D'),5),'k','LineWidth',2);
xlabel('time (min)'); ylabel('Mean FR'); ylim([-1 3]); xlim([0 size(D,1)*Binsize/1e4/60]);
box off;
makepretty

saveas(fh, [directory '/wave_clus/raster_figures/' 'binned_raster_heatmap'], 'png')
saveas(fh, [directory '/wave_clus/raster_figures/' 'binned_raster_heatmap'], 'svg')

close 

%% FIGURE: Matrix of channels cross-correlation
figure('Name','Channels matrix');
imagesc(corr(D))

%% SECTION 3.2: Gamma-FR covariance
% Interpolate firing rate on Gamma's timeline
Mean_FR_on_Gamma = interp1(linspace(0,1,length(D)), movmean(nanmean(D'),5), ...
                   linspace(0,1,length(SmoothGamma)));
Mean_FR_on_Gamma_tsd = tsd(Range(SmoothGamma), Mean_FR_on_Gamma');


Q_sleep = tsd(Range(Restrict(Q,Sleep)),nanzscore(full(Data(Restrict(Q,Sleep)))));
D_Sleep = Data(Q_sleep);

%% FIGURE: Gamma vs Firing Rate
figure('Name','Gamma vs Firing Rate','Color','w');
subplot(2,1,1)
clear R D, R = Range(SmoothGamma,'s'); D = Data(SmoothGamma);
plot(R(1:1e3:end)/60, runmean(D(1:1e3:end),100), 'k');
ylabel('OB gamma power (log scale)'); 
yyaxis right
clear R D, R = Range(Mean_FR_on_Gamma_tsd,'s'); D = Data(Mean_FR_on_Gamma_tsd);
plot(R(1:1e4:end)/60, runmean(D(1:1e4:end),10), 'r');
xlabel('time (min)'); ylabel('FR (zscore)'); xlim([0 240]); makepretty

subplot(2,4,6)
X = Data(Restrict(Mean_FR_on_Gamma_tsd,Wake));
Y = log10(Data(Restrict(SmoothGamma,Wake)));
PlotCorrelations_BM(X(1:1e4:end), Y(1:1e4:end));
xlabel('Mean FR (zscore)'); ylabel('OB gamma power (log scale)');
title('Wake'); axis square; makepretty_BM2; xlim([-1 2.5]); ylim([2.2 2.9]);

subplot(2,4,7)
X = Data(Restrict(Mean_FR_on_Gamma_tsd,Sleep));
Y = log10(Data(Restrict(SmoothGamma,Sleep)));
PlotCorrelations_BM(X(1:1e4:end), Y(1:1e4:end));
xlabel('Mean FR (zscore)'); ylabel('OB gamma power (log scale)');
title('Sleep'); axis square; makepretty_BM2; xlim([-1 2.5]); %ylim([15 3]);

%% FIGURE: FR-Gamma time correlation
figure('Name','FR-Gamma time correlation','Color','w');
[c,lags] = xcorr(X , Y);
plot(lags*(median(diff(Range(Mean_FR_on_Gamma_tsd,'s'))))/60,c)
vline(0,'--r')

%% SECTION 3.3: Correlation Analyses (Gamma, Firing Rates, etc.)
% Firing rate across states: 
Q_Wake = Restrict(Q, Wake); D_Wake = Data(Q_Wake);
Q_NREM = Restrict(Q, SWSEpoch); D_NREM = Data(Q_NREM);
Q_REM = Restrict(Q, REMEpoch); D_REM = Data(Q_REM);
Q_NREM1 = Restrict(Q, (Sleep & Epoch_S2) - REMEpoch); D_NREM1 = Data(Q_NREM1);
Q_NREM2 = Restrict(Q, (Sleep & Epoch_S1) - REMEpoch); D_NREM2 = Data(Q_NREM2);

FiringRate_State(1,:) = nanmean(D);
FiringRate_State(2,:) = nanmean(D_Wake);
FiringRate_State(3,:) = nanmean(D_NREM);
FiringRate_State(4,:) = nanmean(D_REM);
FiringRate_State(5,:) = nanmean(D_NREM1);
FiringRate_State(6,:) = nanmean(D_NREM2);

%% FIGURE: Firing rate across states and channels
figure('Name','Firing Rate Across States','Color','w');
plot(FiringRate_State(1,:), '--k'); hold on
plot(FiringRate_State(2,:), 'b');
plot(FiringRate_State(3,:), 'r');
plot(FiringRate_State(4,:), 'g');
plot(FiringRate_State(5,:),'Color',[1 0.5 0.5]);
plot(FiringRate_State(6,:),'Color',[0.5 0 0]);
xlabel('Neurons'); ylabel('Mean firing rate'); axis square; box off;
makepretty
xlim([1 length(chcl)])

%% FIGURE: Matrix of time correlation across states
figure
subplot(231)
imagesc(corr(D(:,:)')), axis xy, colormap redblue
title('full')
axis square
caxis([-1 1])

subplot(232)
imagesc(corr(D_Wake(:,:)')), axis xy, colormap redblue
title('Wake')
axis square
caxis([-1 1])

subplot(233)
imagesc(corr(D_REM(:,:)')), axis xy, colormap redblue
title('REM')
axis square
caxis([-1 1])

subplot(234)
imagesc(corr(D_NREM(:,:)')), axis xy, colormap redblue
title('NREM')
axis square
caxis([-1 1])

subplot(235)
imagesc(corr(D_NREM1(:,:)')), axis xy, colormap redblue
title('NREM1')
axis square
caxis([-1 1])

subplot(236)
imagesc(corr(D_NREM2(:,:)')), axis xy, colormap redblue
title('NREM2')
axis square
caxis([-1 1])

%% FIGURE: Concatenated time-correlation matrix
figure
imagesc(corr([D_Wake(:,:) ; D_NREM(:,:) ; D_REM(:,:)]')), axis xy, colormap redblue
axis square, caxis([-1 1])
vline([248 997],'-k'), hline([248 997],'-k')

%% FIGURE: FR correlations across states (scatter)
figure('name', 'FR correlations across states')
subplot(331)
PlotCorrelations_BM(FiringRate_State(2,:) , FiringRate_State(3,:))
axis square, makepretty
xlabel('FR Wake'), ylabel('FR NREM')

subplot(332)
PlotCorrelations_BM(FiringRate_State(2,:) , FiringRate_State(4,:))
axis square, makepretty
xlabel('FR Wake'), ylabel('FR NREM')

subplot(333)
PlotCorrelations_BM(FiringRate_State(3,:) , FiringRate_State(4,:))
axis square, makepretty
xlabel('FR NREM'), ylabel('FR REM')


subplot(334)
PlotCorrelations_BM(FiringRate_State(5,:) , FiringRate_State(2,:))
axis square, makepretty
xlabel('FR NREM1'), ylabel('FR Wake')

subplot(335)
PlotCorrelations_BM(FiringRate_State(5,:) , FiringRate_State(4,:))
axis square, makepretty
xlabel('FR NREM1'), ylabel('FR REM')

subplot(336)
PlotCorrelations_BM(FiringRate_State(5,:) , FiringRate_State(6,:))
axis square, makepretty
xlabel('FR NREM1'), ylabel('FR NREM2')

subplot(337)
PlotCorrelations_BM(FiringRate_State(6,:) , FiringRate_State(2,:))
axis square, makepretty
xlabel('FR NREM1'), ylabel('FR Wake')

subplot(338)
PlotCorrelations_BM(FiringRate_State(6,:) , FiringRate_State(4,:))
axis square, makepretty
xlabel('FR NREM1'), ylabel('FR REM')

%% FIGURE: FR correlations across states (matrix)
figure
imagesc(corr(FiringRate_State([2 4:6],:)'))
axis square
xticks([1:5]), yticks([1:5]), xtickangle(45)
xticklabels({'Wake','REM','NREM1','NREM2'}), yticklabels({'Wake','REM','NREM1','NREM2'})
colormap redblue
caxis([-1 1])

%% FIGURE: Box plots for average firing rates in all states
Cols = {[0 0 1],[0 1 0],[1 .5 .5],[.5 0 0]};    
X = [1:4];
Legends = {'Wake','REM','NREM1','NREM2'}; 

figure
MakeSpreadAndBoxPlot3_SB({FiringRate_State(2,:) FiringRate_State(4,:) ...
    FiringRate_State(5,:) FiringRate_State(6,:)},Cols,X,Legends,'showpoints',1,'paired',1)
% set(gca , 'Yscale','log')
ylabel('firing rate (Hz)')
makepretty_BM2

%% FIGURE: Box plots for average firing rates in NREM1 and NREM2
Cols = {[1 .5 .5],[.5 0 0]};
X = [1 2];
Legends = {'NREM1','NREM2'};

figure
MakeSpreadAndBoxPlot3_SB({FiringRate_State(5,:) FiringRate_State(6,:)},Cols,X,Legends,'showpoints',1,'paired',1)
title(['n_{clusters} = ' num2str(length(FiringRate_State))])

%% FIGURE: Histograms
figure
i=1;
subplot(122)
clear d, d=D_NREM2(:,i); d(or(d<-2.5 , d>2.5))=NaN;
hist(d,100)
xlim([-2.5 2.5]);

subplot(121)
clear d, d=D_NREM1(:,i); d(or(d<-2.5 , d>2.5))=NaN;
hist(d,100)
xlim([-2.5 2.5]);

%% SECTION 3.4: Restrict Gamma on spikes across states
SmoothGamma_onSpikes = Restrict(SmoothGamma , Q);
SmoothGamma_onSpikes_Wake = Restrict(SmoothGamma_onSpikes , Wake);
SmoothGamma_onSpikes_NREM = Restrict(SmoothGamma_onSpikes , SWSEpoch);
SmoothGamma_onSpikes_REM = Restrict(SmoothGamma_onSpikes , REMEpoch);
SmoothGamma_onSpikes_NREM1 = Restrict(SmoothGamma_onSpikes , and(Sleep , Epoch_S2)-REMEpoch);
SmoothGamma_onSpikes_NREM2 = Restrict(SmoothGamma_onSpikes , and(Sleep , Epoch_S1)-REMEpoch);

%% FIGURE: Gamma on spikes across states (iterate channels manually)
figure
i = 1;
clf
D1 = log10(Data(SmoothGamma_onSpikes)); D2 = D(:,i); %D3(D2>0) = sqrt(D2(D2>0)); D3(D2<0) = sqrt(-D2(D2<0));
PlotCorrelations_BM(D1 , D2, 'color', [0 0 0] , 'Marker_Size' , 15);
title(['Full, channel #' num2str(metadata(colIndex(i), 1)) ' cluster #' num2str(metadata(colIndex(i), 2))])
i=i+1;

figure
i = 1;
clf
D1 = log10(Data(SmoothGamma_onSpikes_Wake)); D2 = D_Wake(:,i); %D3(D2>0) = sqrt(D2(D2>0)); D3(D2<0) = sqrt(-D2(D2<0));
PlotCorrelations_BM(D1 , D2, 'color', [0 0 1]);
title(['Wake, channel #' num2str(metadata(colIndex(i), 1)) ' cluster #' num2str(metadata(colIndex(i), 2))])
xlabel('OB gamma power (log scale)'), ylabel('FR (zscore)')
i=i+1;

figure
i = 1;
clf
D1 = log10(Data(SmoothGamma_onSpikes_NREM)); D2 = D_NREM(:,i); %D3(D2>0) = sqrt(D2(D2>0)); D3(D2<0) = sqrt(-D2(D2<0));
PlotCorrelations_BM(D1 , D2, 'color', [1 0 0]);
title(['NREM, channel #' num2str(metadata(colIndex(i), 1)) ' cluster #' num2str(metadata(colIndex(i), 2))])
i=i+1;

figure
i = 1;
clf
D1 = log10(Data(SmoothGamma_onSpikes_REM)); D2 = D_REM(:,i); %D3(D2>0) = sqrt(D2(D2>0)); D3(D2<0) = sqrt(-D2(D2<0));
PlotCorrelations_BM(D1 , D2, 'color', [0 1 0]);
title(['REM, channel #' num2str(metadata(colIndex(i), 1)) ' cluster #' num2str(metadata(colIndex(i), 2))])
i=i+1;

figure
i = 1;
clf
D1 = log10(Data(SmoothGamma_onSpikes_NREM1)); D2 = D_NREM1(:,i); %D3(D2>0) = sqrt(D2(D2>0)); D3(D2<0) = sqrt(-D2(D2<0));
PlotCorrelations_BM(D1 , D2, 'color', [1 .5 .5]);
title(['NREM1, channel #' num2str(metadata(colIndex(i), 1)) ' cluster #' num2str(metadata(colIndex(i), 2))])
i=i+1;

figure
i = 1;
clf
D1 = log10(Data(SmoothGamma_onSpikes_NREM2)); D2 = D_NREM2(:,i); %D3(D2>0) = sqrt(D2(D2>0)); D3(D2<0) = sqrt(-D2(D2<0));
PlotCorrelations_BM(D1 , D2, 'color', [.5 0 0]);
title(['NREM2, channel #' num2str(metadata(colIndex(i), 1)) ' cluster #' num2str(metadata(colIndex(i), 2))])
i=i+1;

%% FIGURE: Gamma on spikes across states scatter correlations
figure
for i=1:16
    clf
    subplot(151)
    D1 = log10(Data(SmoothGamma_onSpikes)); D2 = D(:,i); %D3(D2>0) = sqrt(D2(D2>0)); D3(D2<0) = sqrt(-D2(D2<0));
    [R(i) , P(i)] = PlotCorrelations_BM(D1 , D2);
    xlabel('gamma power (log scale)'), ylabel('firing rate (zscore)')
    title('all')
    axis square, box off
    
    subplot(152)
    D1 = log10(Data(SmoothGamma_onSpikes_Wake)); D2 = D_Wake(:,i); D2(D2==0)=NaN;
    [R_Wake(i) , P_Wake(i)] = PlotCorrelations_BM(D1 , D2 , 'color' , [0 0 1]);
    xlabel('gamma power (log scale)'), ylabel('firing rate (zscore)')
    title('Wake')
    axis square, box off
    
    subplot(153)
    D1 = log10(Data(SmoothGamma_onSpikes_NREM1)); D2 = D_NREM1(:,i); D2(D2==0)=NaN;
    [R_NREM1(i) , P_NREM1(i)] = PlotCorrelations_BM(D1 , D2 , 'color' , [1 0 0]);
    xlabel('gamma power (log scale)'), ylabel('firing rate (zscore)')
    title('NREM ')
    axis square, box off
        
    subplot(154)
    D1 = log10(Data(SmoothGamma_onSpikes_NREM2)); D2 = D_NREM2(:,i); D2(D2==0)=NaN;
    [R_NREM2(i) , P_NREM2(i)] = PlotCorrelations_BM(D1 , D2 , 'color' , [1 0 0]);
    xlabel('gamma power (log scale)'), ylabel('firing rate (zscore)')
    title('NREM ')
    axis square, box off
    
    subplot(155)
    D1 = log10(Data(SmoothGamma_onSpikes_REM)); D2 = D_REM(:,i); D2(D2==0)=NaN;
    [R_REM(i) , P_REM(i)] = PlotCorrelations_BM(D1 , D2 , 'color' , [0 1 0]);
    xlabel('gamma power (log scale)'), ylabel('firing rate (zscore)')
    title('REM')
    axis square, box off
end

%% FIGURE: Correlations (r_values)
figure
MakeSpreadAndBoxPlot3_SB({R_Wake(P_Wake<.05) R_NREM1(P_NREM1<.05) R_NREM2(P_NREM2<.05) R_REM(P_REM<.05)},Cols,X,Legends,'showpoints',1,'paired',0)
hline(0,'--r')
ylabel('R values'), ylim([-1 1])


[h,p] = ttest(R_Wake(P_Wake<.05) , zeros(1,12))

%% SECTION 4.1: Build Stimulus Onset States for PSTH
stim_onset = ts(Starttime);
stim_onset_Wake = Restrict(stim_onset, Wake);
stim_onset_Sleep = Restrict(stim_onset, Sleep);
stim_onset_NREM = Restrict(stim_onset, SWSEpoch);
stim_onset_REM = Restrict(stim_onset, REMEpoch);

stim_onset_NREM1 = Restrict(stim_onset, (Sleep & Epoch_S2) - REMEpoch);
stim_onset_NREM2 = Restrict(stim_onset, (Sleep & Epoch_S1) - REMEpoch);

stim_states = {stim_onset_Wake; stim_onset_Sleep; stim_onset_NREM; ...
               stim_onset_REM; stim_onset_NREM1; stim_onset_NREM2};
conditionNames = {'Wake','Sleep','NREM','REM','NREM1','NREM2'};

%% SECTION 4.2: PSTH ANALYSIS (Raw & Z-score) + Population Figures + Controls
fig_vis = 'on';
[rawRates, popRawPSTH, zRates, popZscoredPSTH, zBaselineCorrectedRates, popZscoredPSTH_BC, timeCenters] = PSTH_Analysis_AllNeurons(A, chcl, stim_states, conditionNames, directory, fig_vis);

outFile = fullfile([directory '/wave_clus/'], sprintf('spike_analysis_%s.mat', session));
save(outFile, ...
    'rawRates', 'zRates', 'zBaselineCorrectedRates',...
    'popRawPSTH', 'popZscoredPSTH', 'popZscoredPSTH_BC', ...
    'timeCenters', 'conditionNames', ...
    'metadata', 'colIndex', 'chcl', ...
    '-v7.3');
fprintf('Saved PSTH results to %s\n', outFile);

%% this is just selected relevant blocs from this scrips to run calculate spike_analysis file for all sessions
run_PSTH_analysis_all_sess

%% FIGURE: Simple PSTH with rasters for different states
i = 2; % select the neuron
for j = 1:6
    f7 = figure('Visible', 'off');
    sgtitle(conditionNames{j}, 'FontWeight', 'bold')
    [fh,sq,sweeps] = RasterPETH(A{i}, stim_states{j}, -10000, +50000,'BinSize',240);
    saveas(f7, [directory '/wave_clus/raster_figures/' 'Raster_PSTH_' num2str(i) '_' conditionNames{j}], 'png')
end

%% FIGURE: WakeGamma-dependent PSTH
[rawRatesGamma, popRawPSTHGamma, zRatesGamma, popZscoredPSTHGamma, zBaselineCorrectedRatesGamma, popZscoredPSTH_BCGamma, timeCentersGamma, GammaEpochs, GammaConditionNames] = PSTH_GammaBins_SingleNeuron(A, count, SmoothGamma, Wake, Starttime, directory);

save(outFile, ...
    'rawRatesGamma', 'zRatesGamma', 'zBaselineCorrectedRatesGamma',...
    'popRawPSTHGamma', 'popZscoredPSTHGamma', 'popZscoredPSTH_BCGamma', ...
    'timeCentersGamma', 'GammaEpochs', 'GammaConditionNames', ...
    'metadata', 'colIndex', 'chcl', ...
    '-v7.3', '-append');
fprintf('Saved PSTH results to %s\n', outFile);


%% SECTION 4.3: Take all sessions and plot average figures
MultiSession_Analysis

%% LEGACY: for KlustaKwik preprocessing generate list of .spk files
% name = 'M4_20241206_Shropshire_20241206_fm_TORCs' ;
% % [32 34 35 38 45:55 57 59:72 94 95]
% diary('klusta_input')
% for i = 1:30
%     disp(['KlustaKwik ' name ' ' num2str(i)])
% end
% diary('klusta_input')
% Access KlustaKwik data
% channels = [32 34 35 38 50 51 52 53 54 55 57 60 61 62 63 64 65 66 67 68 69 70 71 72 93 94 45 46 47 48 49 59];  % for 05
% % channels = [31 32 34 35 38 39 45 46 47 49 50 51 52 53 54 55 57 59 60 61 62 68 69 70 82 85 91 92 93 94]; %for 06/09/10/11/12/13/14/21/23 FM TORCs
% % ACx1 = [31 32 34 35 38 39 45 46 47 49 50 51 52 53 54 55 57 59 60 61 62]
% % ACx2 = [68 69 70 82 85 91 92 93 94]
% 
% cd([directory '/SpikeSorting'])
% 
% SetCurrentSession(name)
% % SetCurrentSession('same')
% spikes = GetSpikes('all','output','full');

% For KlustaKwik
% for ii = 1:size(channels)
%     figure
% %     [fh,sq,sweeps] = RasterPETH(ts(spikes(spikes(:,2)==ii,1)*1e4), ts(Starttime*2/3), -10000, +15000,'BinSize',100);
%     [fh,sq,sweeps] = RasterPETH(ts(spikes(spikes(:,2)==ii,1)*1e4), ts(Starttime), -10000, +50000,'BinSize',240);
%     title(['channel #' num2str(channels(ii))])
%     Resp(ii,:) = Data(sq);
% end
%{
% name = 'M4_20241206_Shropshire_20241206_fm_TORCs' ;
% diary('klusta_input')
% for i = 1:30
%     disp(['KlustaKwik ' name ' ' num2str(i)])
% end
% diary('klusta_input')
%}


end