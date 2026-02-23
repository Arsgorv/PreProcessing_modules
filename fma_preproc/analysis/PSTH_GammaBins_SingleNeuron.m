function [rawRates, popRawPSTH, zRates, popZscoredPSTH, zBaselineCorrectedRates, popZscoredPSTH_BC, timeCenters, binEpochs, conditionNames] = PSTH_GammaBins_SingleNeuron(A, chcl, SmoothGamma, Wake, Starttime, directory, fig_vis)
% PSTH_GammaBins_SingleNeuron
%
% 1) Restrict gamma to Wake: Gamma_Wake
% 2) Bin the gamma power into 10 equal bins (lowest->highest)
% 3) For each bin, restrict stimulus onset times to that bin
% 4) Compute PSTH for each bin
% 5) Plot 10 PSTHs on the same figure, color-coded from blue (low gamma) to red (high gamma)
%
% INPUTS:
%   A         = cell array of ts objects for your spikes
%   iNeuron   = index of the neuron in A (e.g. 2)
%   SmoothGamma = tsd with gamma power (time vs power)
%   Wake      = intervalSet (the Wake epochs)
%   Starttime = stimulus onsets (ts object), or numeric times
%
% By: A. Goriachenkov 02/2025
% ----------------------------------------------------------

%
% Gamma_Wake = Restrict(SmoothGamma, Wake);
% [Y,X]=hist(log(Data(Gamma_Wake)),1000);
% plot(X,Y)
% 
% 
% Data(Gamma_Wake)
% 
% figure
% [Y,X] = hist(log10(Data(Gamma_Wake)),1000);
% a = plot(X , runmean(Y,3));
% % a.FaceColor=[.8 .8 .8]; a.LineWidth=3; a.EdgeColor=[0 0 0];
% box off

%%
conditionNames = {'Lowest Gamma', 'Low Gamma', 'High Gamma', 'Highest gamma'};
allCondsIdx = 1:length(conditionNames);

% PSTH window
t_pre  = -1;
t_post =  5;
binSize= 0.05; % 50 ms bins

timeEdges   = t_pre : binSize : t_post;
timeCenters = (timeEdges(1:end-1) + timeEdges(2:end))/2;

baselineIdx = find(timeCenters >= -1 & timeCenters <= 0);
peakIdx = find(timeCenters >= 0 & timeCenters <= 0.1);
ongoingIdx = find(timeCenters >= 0.075 & timeCenters <= 3);
stimIdx = find(timeCenters >= 0 & timeCenters <= 3);

% Stim highlight
baselineWindow = [timeCenters(baselineIdx(1)) timeCenters(baselineIdx(end))];
peakWindow = [timeCenters(peakIdx(1)) timeCenters(peakIdx(end))];
ongoingWindow = [timeCenters(ongoingIdx(1)) timeCenters(ongoingIdx(end))];
stimWindow = [timeCenters(stimIdx(1)) timeCenters(stimIdx(end))];

% Colors for each condition
cmap = parula(length(allCondsIdx));  % or your choice: jet(10), hot(10), etc.
myColors = {cmap(1, :), cmap(2, :), cmap(3, :), cmap(4, :)};

yminZ=-2; ymaxZ=4; 

%% 1) Initialize population arrays
nNeurons = length(A);
nConds = length(allCondsIdx);
nBins = length(timeCenters);

popRawPSTH  = nan(nNeurons, nBins, nConds);
popZscoredPSTH = nan(nNeurons, nBins, nConds);
popZscoredPSTH_BC = nan(nNeurons, nBins, nConds);
popAvgRate = nan(nNeurons, nConds);
popAvgZ = nan(nNeurons, nConds);

%% 2) Define binsize
stim_onset = ts(Starttime); 
Wake_stim = Restrict(stim_onset, Wake);

%% 3) Restrict gamma to Wake
Gamma_Wake = Restrict(SmoothGamma, Wake);

% For each stim onset, find the exact gamma value
pStim = interp1(Range(Gamma_Wake), Data(Gamma_Wake), Range(Wake_stim), 'linear');

edges = quantile(pStim, [0 0.25 0.5 0.75 1]);

% binsize = round(length(Range(Wake_stim))/4); %ensure equal number of trials in each bin

%% 4) We'll create intervalSets for each bin. The logic:
%   For bin i, we keep times tGamma where pGamma is in [edges(i), edges(i+1))
GammaBins = 4;
binEpochs = cell(GammaBins,1);

% Bin 1
Epoch_1 = thresholdIntervals(Gamma_Wake , edges(1), 'Direction','Above');
Epoch_2 = thresholdIntervals(Gamma_Wake , edges(2), 'Direction','Below');
Epoch_G_low = and(Epoch_1, Epoch_2);
% Gamma_low = Restrict(Gamma_Wake, Epoch_G_low);

% Bin 2
clear Epoch_1 Epoch_2 Epoch
Epoch_1 = thresholdIntervals(Gamma_Wake , edges(2), 'Direction','Above');
Epoch_2 = thresholdIntervals(Gamma_Wake , edges(3), 'Direction','Below');
Epoch_mid_low = and(Epoch_1, Epoch_2);
% Gamma_mid_low = Restrict(Gamma_Wake, Epoch_mid_low);

% Bin 3
clear Epoch_1 Epoch_2 Epoch
Epoch_1 = thresholdIntervals(Gamma_Wake , edges(3), 'Direction','Above');
Epoch_2 = thresholdIntervals(Gamma_Wake , edges(4), 'Direction','Below');
Epoch_mid_high = and(Epoch_1, Epoch_2);
% Gamma_mid_high = Restrict(Gamma_Wake, Epoch_mid_high);

% Bin 3
clear Epoch_1 Epoch_2 Epoch
Epoch_1 = thresholdIntervals(Gamma_Wake , edges(4), 'Direction','Above');
Epoch_2 = thresholdIntervals(Gamma_Wake , edges(5), 'Direction','Below');
Epoch_high = and(Epoch_1, Epoch_2);
% Gamma_high = Restrict(Gamma_Wake, Epoch_high);

binEpochs = {Epoch_G_low; Epoch_mid_low; Epoch_mid_high; Epoch_high};

%% 5) For each bin, restrict the stimulus onset times to that bin
% We'll store in cell array
stim_states = cell(GammaBins,1);

for i = 1:GammaBins
    stim_states{i} = Restrict(Wake_stim, binEpochs{i});
end

%% 6) Plot Figures

fprintf('plotting PSTHs\n');
trialLabels = [];
for cn = 1:length(stim_states)
    trialLabels = [trialLabels cn*ones(1,length(Range(stim_states{cn})))];
end
for iNeuron=1:nNeurons
    clear f1
    f1 = figure('Color','w','Name',sprintf('Neuron %d PSTH', iNeuron),...
        'Units','normalized','Position',[0.1 0.3 0.7 0.4],'Visible', fig_vis);

    rawRates = cell(nConds,1);
    zRates = cell(nConds,1);

    axZ = subplot(1,2,1); hold on; title('Z-scored PSTH'); 
    xlabel('Time (s)'); ylabel('Z-score');
    axZ_BC = subplot(1,2,2); hold on; title('Z-scored Baseline-corrected PSTH'); 
    xlabel('Time (s)'); ylabel('Z-score');    

    % Gray patch
    patch(axZ_BC,[stimWindow(1) stimWindow(2) stimWindow(2) stimWindow(1)],...
        [yminZ yminZ ymaxZ ymaxZ],[0.8 0.8 0.8],'EdgeColor','none',...
        'FaceAlpha',0.5,'HandleVisibility','off');
    patch(axZ,[stimWindow(1) stimWindow(2) stimWindow(2) stimWindow(1)],...
        [yminZ yminZ ymaxZ ymaxZ],[0.8 0.8 0.8],'EdgeColor','none',...
        'FaceAlpha',0.5,'HandleVisibility','off');

    lineHandlesZ = gobjects(nConds,1);
    lineHandlesRaw = gobjects(nConds,1);

    tmp_raw_acrossC = [];
    for cIdx=1:nConds
        c = allCondsIdx(cIdx);
        onsetTimes = Range(stim_states{c},'s');
        nOnsets = length(onsetTimes);

        % Trial-by-trial PSTH
        tmp_raw = zeros(nOnsets, nBins);
        for iOnset=1:nOnsets
            sTimes = Range(A{iNeuron},'s');
            relTimes = sTimes - onsetTimes(iOnset);
            inWindow = (relTimes>=t_pre & relTimes<=t_post);
            theseSpks = relTimes(inWindow);
            counts = histcounts(theseSpks,timeEdges);
            tmp_raw(iOnset,:) = counts / binSize; 
        end
        tmp_raw_acrossC = [tmp_raw_acrossC ; tmp_raw];
    end
    % [neuron, trial, timebin]
    % tmp_raw_all
    % compute mean for each stage across stages, and std aross stages, by
    % stratifying across conditions
    for cn = 1:max(trialLabels) % min. nb of trial per condition
        trialPerC(cn) = sum(trialLabels==cn);
    end
    minTrialPerC = min(trialPerC);
    clear pickedupIdx;
    pickedupIdx.all = [];
    for cn = 1:max(trialLabels)
        idx = randperm(trialPerC(cn),minTrialPerC);
        tl = find(trialLabels==cn);
        idx = tl(idx);
        pickedupIdx.all = [pickedupIdx.all idx];
        pickedupIdx.C{cn} = idx;
    end
    
    idxBL = (timeCenters>=-1 & timeCenters<0);
    baseline = tmp_raw_acrossC(pickedupIdx.all,idxBL);
    meanBL = mean(baseline,1:2); %mean across a subset of trials (same trial nb in each stage)
    stdBL = std(mean(baseline,2),1); %std across a subset of trials (same trial nb in each stage)
    for cn = 1:max(trialLabels)
        meanBLperStage(cn) = mean(tmp_raw_acrossC(trialLabels==cn,idxBL),1:2);
        % mean across all trials of each stage
    end

    for cIdx=1:nConds       
        tmp_raw = tmp_raw_acrossC(trialLabels==cIdx,:);
        tmp_z = (tmp_raw - meanBL) ./ stdBL;
        tmp_z_BC = (tmp_raw - meanBLperStage(cIdx)) ./ stdBL;

        rawRates{cIdx} = tmp_raw;
        zRates{cIdx} = tmp_z;
        zBaselineCorrectedRates{cIdx} = tmp_z_BC;
        popRawPSTH(iNeuron,:,cIdx) = mean(tmp_raw,1);
        popZscoredPSTH(iNeuron,:,cIdx) = mean(tmp_z,1);
        popZscoredPSTH_BC(iNeuron,:,cIdx) = mean(tmp_z_BC,1);

        % Plot: z-scored baseline corrected
        meanZ_BC = mean(tmp_z_BC,1);
        semZ_BC = std(tmp_z_BC,0,1)/sqrt(size(tmp_z_BC,1));
        ciUpZ_BC = meanZ_BC + 1.96*semZ_BC; 
        ciLoZ_BC = meanZ_BC - 1.96*semZ_BC;

        subplot(axZ_BC);
        fill([timeCenters fliplr(timeCenters)], ...
             [ciUpZ_BC fliplr(ciLoZ_BC)], myColors{cIdx}, ...
             'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
        lineHandlesZ(cIdx) = plot(timeCenters,meanZ_BC,'Color',myColors{cIdx}, ...
            'LineWidth',2,'DisplayName',conditionNames{cIdx});

        % Plot: z-scored
        meanZ = mean(tmp_z,1);
        semZ = std(tmp_z,0,1)/sqrt(size(tmp_z,1));
        ciUpZ = meanZ + 1.96*semZ; 
        ciLoZ = meanZ - 1.96*semZ;

        subplot(axZ);
        fill([timeCenters fliplr(timeCenters)], ...
             [ciUpZ fliplr(ciLoZ)], myColors{cIdx}, ...
             'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
        lineHandlesRaw(cIdx) = plot(timeCenters,meanZ,'Color',myColors{cIdx}, ...
            'LineWidth',2,'DisplayName',conditionNames{cIdx});
    end

    legend(axZ_BC, lineHandlesZ, {conditionNames{allCondsIdx}},'Location','best');
    legend(axZ,lineHandlesRaw,{conditionNames{allCondsIdx}},'Location','best');
    xlim(axZ_BC,[t_pre t_post]); ylim(axZ_BC,[yminZ ymaxZ]); box(axZ_BC,'off');
    xlim(axZ,[t_pre t_post]); ylim(axZ,[yminZ ymaxZ]); box(axZ,'off');

    saveas(f1, [directory '/wave_clus/raster_figures/' 'PSTHs_gamma_Wake_#' num2str(chcl{iNeuron}(1)) '_' num2str(chcl{iNeuron}(2)) '_bs_' num2str(binSize)], 'png')
end

%% 3.1) Population Figures (z-scored PSTH, z-scored BC PSTH)
fprintf('plotting population PSTH\n');
% (A) Population PSTH
% close all
clear f2
f2 = figure('Color','w','Name','Population PSTH','Units','normalized','Position',[0.1 0.3 0.7 0.4],'Visible', fig_vis);
% Z-scored
ax1=subplot(1,2,1); hold on; title('Population-Avg Z-scored PSTH'); xlabel('Time (s)'); ylabel('Z-score');
patch(ax1,[stimWindow(1) stimWindow(2) stimWindow(2) stimWindow(1)], [yminZ yminZ ymaxZ ymaxZ], ...
    [0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.5,'HandleVisibility','off');
lineHandlesPopRaw = gobjects(nConds,1);
for cIdx=1:nConds
    c = allCondsIdx(cIdx);
    rawAllNeurons = squeeze(popZscoredPSTH(:,:,cIdx));
    meanRaw_pop = mean(rawAllNeurons,1);
    semRaw_pop = std(rawAllNeurons,0,1)/sqrt(size(rawAllNeurons, 1));
    ciUp = meanRaw_pop + 1.96*semRaw_pop;
    ciLo = meanRaw_pop - 1.96*semRaw_pop;
    fill([timeCenters fliplr(timeCenters)], [ciUp fliplr(ciLo)], ...
         myColors{c}, 'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
    lineHandlesPopRaw(cIdx) = plot(timeCenters, meanRaw_pop, 'Color', myColors{c},...
        'LineWidth',2,'DisplayName',conditionNames{c});
end
legend(ax1, lineHandlesPopRaw,{conditionNames{allCondsIdx}},'Location','best');
xlim(ax1,[t_pre t_post]); ylim(ax1,[-1 2.5]); box(ax1,'off');

% Z-scored BC
ax2=subplot(1,2,2); hold on; title('Population-Avg Z-scored baseline-corrected PSTH'); xlabel('Time (s)'); ylabel('Z-score');
patch(ax2,[stimWindow(1) stimWindow(2) stimWindow(2) stimWindow(1)], [yminZ yminZ ymaxZ ymaxZ], ...
    [0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.5,'HandleVisibility','off');
lineHandlesPopZ = gobjects(nConds,1);
for cIdx=1:nConds
    c = allCondsIdx(cIdx);
    zAllNeurons = squeeze(popZscoredPSTH_BC(:,:,cIdx));
    meanZ_pop = mean(zAllNeurons,1);
    semZ_pop = std(zAllNeurons,0,1)/sqrt(size(zAllNeurons, 1));
    ciUpz = meanZ_pop + 1.96*semZ_pop;
    ciLoz = meanZ_pop - 1.96*semZ_pop;
    fill([timeCenters fliplr(timeCenters)], [ciUpz fliplr(ciLoz)], ...
         myColors{c}, 'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
    lineHandlesPopZ(cIdx) = plot(timeCenters, meanZ_pop, 'Color', myColors{c},...
        'LineWidth',2,'DisplayName',conditionNames{c});
end
legend(ax2, lineHandlesPopZ,{conditionNames{allCondsIdx}},'Location','best');
xlim(ax2,[t_pre t_post]); ylim(ax2,[yminZ ymaxZ]); box(ax2,'off');

saveas(f2, [directory '/wave_clus/raster_figures/' 'population_PSTHs_gamma_Wake_bs_' num2str(binSize)], 'png')

%% 3.2) Population Figures (Box plots)
clear f3

f3 = figure('Color','w','Name','Population Box','Units','normalized','Position',[0.1 0.3 0.7 0.4],'Visible', fig_vis);

ax1=subplot(1,3,1); hold on; title(['Avg baseline rate across neurons, n = ' num2str(nNeurons)]); ylabel('Z-score');
ax2=subplot(1,3,2); hold on; title(['Avg peak-response across neurons, n = ' num2str(nNeurons)]); ylabel('Z-score');
ax3=subplot(1,3,3); hold on; title(['Avg ongoing response across neurons, n = ' num2str(nNeurons)]); ylabel('Z-score');

myData=cell(3,1); myGroups = cell(3,1);

popAvgRate_baseline = squeeze(mean(popZscoredPSTH(:, find(timeCenters >= baselineWindow(1) & timeCenters <= baselineWindow(2)), :), 2));
popAvgRate_peak = squeeze(mean(popZscoredPSTH_BC(:, find(timeCenters >= peakWindow(1) & timeCenters <= peakWindow(2)), :), 2));
popAvgRate_ongoing = squeeze(mean(popZscoredPSTH_BC(:, find(timeCenters >= ongoingWindow(1) & timeCenters <= ongoingWindow(2)), :), 2));

popAvgRate = {popAvgRate_baseline; popAvgRate_peak; popAvgRate_ongoing};

axs = {ax1, ax2, ax3};
for it = 1:3
    for cIdx=1:nConds
        theseVals{it} = popAvgRate{it}(:,cIdx);
        myData{it} = [myData{it}; theseVals{it}];
        myGroups{it} = [myGroups{it}; cIdx*ones(size(theseVals{it}))];
    end
    boxplot(axs{it}, myData{it}, myGroups{it}, 'Labels', conditionNames(allCondsIdx));
    xlim(axs{it},[0.5 nConds+0.5]); box(axs{it},'off');
end
xtickangle(ax1, 45);
xtickangle(ax2, 45);
xtickangle(ax3, 45);
% [p,~,stats] = anova1(myData,myGroups,'off');
% fprintf('Population-level ANOVA p=%.4g\n', p);
% posthoc=multcompare(stats,'ctype','tukey-kramer'); disp(posthoc);

saveas(f3, [directory '/wave_clus/raster_figures/' 'population_box_plots_gamma_Wake_bs_' num2str(binSize)], 'png')

%% Baptiste's Box plots
clear f3_1

f3_1 = figure('Color','w','Name','Population Box','Units','normalized','Position',[0.1 0.3 0.7 0.4],'Visible', fig_vis);

ax1_1=subplot(1,3,1); hold on; title(['Avg baseline rate across neurons, n = ' num2str(nNeurons)]); ylabel('Z-score');
ax2_1=subplot(1,3,2); hold on; title(['Avg peak-response across neurons, n = ' num2str(nNeurons)]); ylabel('Z-score');
ax3_1=subplot(1,3,3); hold on; title(['Avg ongoing response across neurons, n = ' num2str(nNeurons)]); ylabel('Z-score');

for i = 1:3
    subplot(1,3,i)
    MakeSpreadAndBoxPlot3_SB(popAvgRate{i}, myColors, [1:4], conditionNames(allCondsIdx),'showpoints',1,'paired',1)
end
saveas(f3_1, [directory '/wave_clus/raster_figures/' 'population_box_plots_gamma_Wake_BM_bs_' num2str(binSize)], 'png')

end
