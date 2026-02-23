function [rawRates, popRawPSTH, zRates, popZscoredPSTH, zBaselineCorrectedRates, popZscoredPSTH_BC, timeCenters] = PSTH_Analysis_AllNeurons(A, chcl, stim_states, conditionNames, directory, fig_vis)
% PSTH_Analysis_AllNeurons
%
% OUTPUTS:
%   popRawPSTH      (nNeurons x nBins x nConds)
%   popZscoredPSTH  (nNeurons x nBins x nConds)
%   popAvgRate      (nNeurons x nConds)
%   popAvgZ         (nNeurons x nConds)
%   timeCenters     (1 x nBins)
%
% Builds raw & z-scored PSTHs for all neurons in A, across all conditions 
% in stim_states. Then:
%   - Plots per-neuron PSTH (2-subplot: raw + z-score)
%   - Population PSTH (raw + z)
%   - Violin plots of average firing
%   - Heatmap of z-scored PSTH
%   - Control analyses: sub-sampling, random onsets, spike jitter
%
% USAGE:
%   PSTH_Analysis_AllNeurons(A, stim_states, conditionNames, directory)

fprintf('== Running PSTH_Analysis_AllNeurons ==\n');

if length(conditionNames) ~= length(stim_states)
    error('conditionNames and stim_states must have the same length!');
end

%% 0) Config
allCondsIdx = 1:length(conditionNames);

% PSTH window
t_pre  = -1;
t_post =  5;
binSize= 0.05; % 50 ms bins

% timeEdges = (t_pre - binSize/2) : binSize : (t_post + binSize/2);
% 
% timeCenters = timeEdges(1:end-1) + binSize/2;

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
myColors = {
    [0 0 1];     % Wake
    [.7 .7 .7];     % Sleep
    [1 0 0];   % NREM
    [0 1 0];     % REM
    [1 0.5 0.5]; % NREM1
    [0.5 0 0];   % NREM2
};

yminZ=-1; ymaxZ=4; 

%% 1) Initialize population arrays
nNeurons = length(A);
nConds = length(allCondsIdx);
nBins = length(timeCenters);

popRawPSTH  = nan(nNeurons, nBins, nConds);
popZscoredPSTH = nan(nNeurons, nBins, nConds);
popZscoredPSTH_BC = nan(nNeurons, nBins, nConds);
popAvgRate = nan(nNeurons, nConds);
popAvgZ = nan(nNeurons, nConds);

%% 2) Loop over Neurons & Build/Plot PSTHs
fprintf('plotting PSTHs\n');
% tmp_raw_all = [];
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
    xlim(axZ_BC,[t_pre t_post]); ylim(axZ_BC,[-1 3]); box(axZ_BC,'off');
    xlim(axZ,[t_pre t_post]); ylim(axZ,[-1 3]); box(axZ,'off');

    saveas(f1, [directory '/wave_clus/raster_figures/' 'PSTHs_#' num2str(chcl{iNeuron}(1)) '_' num2str(chcl{iNeuron}(2)) '_bs_' num2str(binSize)], 'png')
end

%% 3.1) Population Figures (z-scored PSTH, z-scored BC PSTH)
fprintf('plotting population PSTH\n');
% (A) Population PSTH
close all
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

saveas(f2, [directory '/wave_clus/raster_figures/' 'population_PSTHs_bs_' num2str(binSize)], 'png')

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
% [p,~,stats] = anova1(myData,myGroups,'off');
% fprintf('Population-level ANOVA p=%.4g\n', p);
% posthoc=multcompare(stats,'ctype','tukey-kramer'); disp(posthoc);

saveas(f3, [directory '/wave_clus/raster_figures/' 'population_box_plots_bs_' num2str(binSize)], 'png')

%% Baptiste's Box plots
clear f3_1

f3_1 = figure('Color','w','Name','Population Box','Units','normalized','Position',[0.1 0.3 0.7 0.4],'Visible', fig_vis);

ax1_1=subplot(1,3,1); hold on; title(['Avg baseline rate across neurons, n = ' num2str(nNeurons)]); ylabel('Z-score');
ax2_1=subplot(1,3,2); hold on; title(['Avg peak-response across neurons, n = ' num2str(nNeurons)]); ylabel('Z-score');
ax3_1=subplot(1,3,3); hold on; title(['Avg ongoing response across neurons, n = ' num2str(nNeurons)]); ylabel('Z-score');

for i = 1:3
    subplot(1,3,i)
    MakeSpreadAndBoxPlot3_SB(popAvgRate{i}, myColors, [1:6], conditionNames(allCondsIdx),'showpoints',0,'paired',0)
end
saveas(f3_1, [directory '/wave_clus/raster_figures/' 'population_box_plots_BM_bs_' num2str(binSize)], 'png')

%% 3.4) Heatmaps (Z-scored)
fprintf('plotting heat maps\n');
% close all
clear f4
f4 = figure('Color','w','Name','PSTH Heatmaps','Units','normalized','Position',[0.2 0.1 0.7 0.7],'Visible', fig_vis);
[nNeurons, ~, nConds] = size(popZscoredPSTH_BC);
for c=1:nConds
    subplot(2,3,c);
    dataMatrix = squeeze(popZscoredPSTH_BC(:,:,c));  % (nNeurons x nBins)
    imagesc(timeCenters, 1:nNeurons, dataMatrix);
    set(gca,'YDir','normal'); colormap(gca,'viridis');
    caxis([-0.3 0.7]); colorbar;
    title(conditionNames{c}); xlabel('Time (s)'); ylabel('Neuron #');
end
sgtitle('Heatmaps of Z-scored BC PSTH across Neurons');

saveas(f4, [directory '/wave_clus/raster_figures/' 'heat_maps_PSTH_' num2str(binSize)], 'png')


%% NOT USED 3.3) Violin Plot
% fprintf('plotting violin plot\n');
% close all
% clear f4
% f4 = figure('Color','w','Name','Violin Plots: Per-Neuron Mean Rate','Visible', fig_vis);
% hold on; title('Violin Plot of [0-3s] Firing Rate across Neurons');
% ylabel('Mean Firing Rate [spikes/s]');
% [nN, nC] = size(popAvgRate);
% allData=[]; allGroups=[];
% for c=1:nC
%     vals = popAvgRate(:,c);
%     allData=[allData; vals]; 
%     allGroups=[allGroups; c*ones(nN,1)];
% end
% violinplot(allData, allGroups,'ShowMean', true); %,'ViolinColor', cell2mat(myColors), 'Labels', cell2mat(conditionNames), 
% box off; hold off;
% 
% saveas(f4, [directory '/wave_clus/raster_figures/' 'violin_FR_across_neurons'], 'png')


%% NOT USED (rewrite ComputePSTH.m: CONTROL ANALYSES (Sub-sampling, random onsets, spike jitter)
% nNeuron=2; % pick an example neuron
% 
% % i) Sub-Sampling
% fprintf('plotting control: subsampling\n');
% close all
% clear f5
% f5 = figure('Name','Sub-sampled PSTH','Color','w','Units','normalized','Position',[0.1 0.3 0.7 0.4],'Visible', fig_vis); 
% sgtitle(['Control: Sub-sampled PSTH for channel#' num2str(nNeuron)], 'FontWeight', 'bold', 'FontSize', 15)
% ax4 = subplot(1,2,2);hold on;
% title(ax4, ['Z-scored PSTH']);
% xlabel(ax4, 'Time (s)'); ylabel(ax4, 'Z-score'); legend(ax4, 'Location','best'); box off;
% ax5 = subplot(1,2,1);hold on;
% xlabel(ax5, 'Time (s)'); ylabel(ax5, 'Spikes/s'); legend(ax5, 'Location','best'); box off;
% title(ax5, ['Raw PSTH']);
% spikeTimes=Range(A{nNeuron},'s');
% trialCounts=zeros(1,nConds);
% for c=1:nConds
%     trialCounts(c)=length(Range(stim_states{c},'s'));
% end
% Nmin=min(trialCounts);
% 
% for c=1:nConds
%     onAll = Range(stim_states{c},'s');
%     idx = randperm(trialCounts(c), Nmin);
%     subsetOnsets = onAll(idx);
%     [mP, sP, ~, mPz, sPz, ~] = ComputePSTH(spikeTimes, subsetOnsets, timeEdges);
%     ciUp = mP + 1.96*sP; ciLo = mP - 1.96*sP;
%     ciUp_z = mPz + 1.96*sPz; ciLo_z = mPz - 1.96*sPz;
%     
% 
%     fill(ax4, [timeCenters fliplr(timeCenters)], [ciUp fliplr(ciLo)], ...
%          myColors{c}, 'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
%     plot(ax4, timeCenters,mP,'Color',myColors{c},'LineWidth',2,'DisplayName',conditionNames{c});
% 
%     fill(ax5, [timeCenters fliplr(timeCenters)], [ciUp_z fliplr(ciLo_z)], ...
%          myColors{c}, 'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
%     plot(ax5, timeCenters,mPz,'Color',myColors{c},'LineWidth',2,'DisplayName',conditionNames{c});
% end
% 
% saveas(f5, [directory '/wave_clus/raster_figures/' 'control_subsampling_PSTH'], 'png')
% 
% 
% % ii) Random Onsets
% fprintf('plotting control: random onsets\n');
% close all
% clear f6
% f6 = figure('Name','Shuffle Control','Color','w','Units','normalized','Position',[0.1 0.3 0.7 0.4],'Visible', 'off');
% sgtitle('Shuffle Control: True PSTH vs Random Onset PSTH', 'FontWeight', 'bold', 'FontSize', 15);
% 
% ax(1) = subplot(1,2,1); hold on; title('True PSTH');
% ax(2) = subplot(1,2,2); hold on; title('Random-Onset PSTH');
% c=1; % e.g. Wake
% 
% trueOnsets=Range(stim_states{c},'s');
% [mean_true,sem_true]=ComputePSTH(spikeTimes,trueOnsets,timeEdges);
% tc=timeCenters;
% subplot(1,2,1);
% fill([tc fliplr(tc)], [mean_true+1.96*sem_true fliplr(mean_true-1.96*sem_true)], ...
%      'b','FaceAlpha',0.2,'EdgeColor','none');
% plot(tc,mean_true,'b','LineWidth',2); xlabel('Time(s)'); ylabel('Spikes/s');
% 
% nOnsets=length(trueOnsets);
% % load([directory '/SleepScoring_OBGamma.mat'], 'Epoch')
% % temp=End(Epoch)/1e4;
% temp=trueOnsets(end);
% randOnsets=sort(rand(nOnsets,1)*temp(end));
% [mean_rand,sem_rand]=ComputePSTH(spikeTimes,randOnsets,timeEdges);
% subplot(1,2,2);
% fill([tc fliplr(tc)], [mean_rand+1.96*sem_rand fliplr(mean_rand-1.96*sem_rand)], ...
%      'r','FaceAlpha',0.2,'EdgeColor','none');
% plot(tc,mean_rand,'r','LineWidth',2); xlabel('Time(s)'); ylabel('Spikes/s');
% linkaxes(ax, 'y')
% saveas(f6, [directory '/wave_clus/raster_figures/' 'control_random_onsets_PSTH'], 'png')
% 
% % iii) Spike Jitter
% fprintf('plotting control: spike jitter\n');
% 
% J=0.5; % ±500 ms
% close all
% clear f7
% f7 = figure('Name','Spike Jitter','Color','w','Visible', 'off'); hold on;
% [mean_orig,sem_orig] = ComputePSTH(spikeTimes,trueOnsets,timeEdges);
% jitterVec = (2*J)*rand(size(spikeTimes)) - J;
% jitSpikes = spikeTimes + jitterVec;
% [mean_jit,sem_jit] = ComputePSTH(jitSpikes,trueOnsets,timeEdges);
% fill([tc fliplr(tc)], [mean_orig+1.96*sem_orig fliplr(mean_orig-1.96*sem_orig)], ...
%      'b','FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
% plot(tc,mean_orig,'b','LineWidth',2,'DisplayName','Original');
% fill([tc fliplr(tc)], [mean_jit+1.96*sem_jit fliplr(mean_jit-1.96*sem_jit)], ...
%      'r','FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
% plot(tc,mean_jit,'r','LineWidth',2,'DisplayName','Jittered');
% xlabel('Time(s)'); ylabel('Spikes/s'); legend('Location','best'); box off;
% title(sprintf('Spike Jitter (Neuron %d, Cond=%s)',nNeuron,conditionNames{c}));
% saveas(f7, [directory '/wave_clus/raster_figures/' 'control_spike_jitter'], 'png')

end