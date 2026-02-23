function run_zeta_test_AG(directory, sessName)
% This function runs zeta test for neruons of a selected session
disp(['running session ' sessName])
dblUseMaxDur = 3.5; % Ignore spikes beyond 3.5s after onset
subsampling_value = 0.5; % 0.5 = 50% of all trials
intResampNum = 1000; % number of resampling

% load spikes
[spikes, metadata] = load_wave_clus(directory);

% load stims
lfpFile = fullfile(directory, 'LFPData', 'LFP111.mat');
if exist(lfpFile,'file')
    load(lfpFile,'LFP');
else
    error('LFP file not found: %s', lfpFile);
end

% Extract stimulus onsets:
StartSound = thresholdIntervals(LFP, 4000, 'Direction','Above');
Starttime  = Start(StartSound); 
% subsample trials
trial_onsets = (Starttime/1e4);
idx = sort(randperm(length(trial_onsets), round(length(trial_onsets)*subsampling_value)))';
vecStimulusStartTimes = trial_onsets(idx);
vecStimulusStopTimes = vecStimulusStartTimes+3;
matEventTimes = cat(2,vecStimulusStartTimes,vecStimulusStopTimes); % put stimulus start and stop times together into a [T x 2] matrix

%% loop over neurons
sigs = nan(1,size(spikes, 2));

for s_count = 1:size(spikes, 2)
    vecSpikeTimes1 = spikes(:, s_count)/1e3; %in seconds
    [dblZetaP] = zetatest(vecSpikeTimes1, matEventTimes, dblUseMaxDur, intResampNum);
    sigs(s_count) = dblZetaP;
    i = metadata(s_count,:);
    fprintf('Channel #%d, Cluster #%d, p-value = %d\n', i(1), i(2), sigs(s_count))
end

signif_idx = find(sigs <= 0.05);
list_signif_trials_1000 = [metadata(signif_idx, :), sigs(signif_idx)'];
list_all_trials_1000 = [metadata, sigs'];
non_signif_idx = find(sigs > 0.05);

%% plot fig
fh = figure('Visible', 'off');
subplot(1,2,1)
bar(signif_idx, sigs(signif_idx), 'FaceColor', 'r')
hold on
bar(non_signif_idx, sigs(non_signif_idx),  'FaceColor', 'k')
xticks(1:size(spikes, 2))

sgtitle(['Session ' sessName])
xlabel('Channel Number')
ylabel('ZETA P val')

subplot(1,2,2)
bar(signif_idx, sigs(signif_idx), 'FaceColor', 'r')
hold on
bar(non_signif_idx, sigs(non_signif_idx),  'FaceColor', 'k')
xticks(1:size(spikes, 2))

xlabel('Channel Number')
ylabel('ZETA P val')
ylim([0 0.05])

%% save stuff
saveas(fh, [directory '/wave_clus/raster_figures/zeta_test_1000'], 'png');

% save([directory '/wave_clus/spike_analysis_' sessName], 'list_all_trials_50', 'list_signif_trials_50', '-append')
save([directory '/wave_clus/spike_analysis_' sessName], 'list_all_trials_1000', 'list_signif_trials_1000')

close(fh);

%% run the ZETA-test with specified parameters
% %however, we can also specify the parameters ourselves
% dblUseMaxDur = min(diff(vecStimulusStartTimes)); %minimum of trial-to-trial durations
% % intResampNum = 50; %~50 random resamplings should give us a good enough idea if this cell is responsive, but if it's close to 0.05, we should increase this #. Generally speaking, more is better, so let's put 100 here.
% intPlot = 3;%what do we want to plot?(0=nothing, 1=inst. rate only, 2=traces only, 3=raster plot as well, 4=adds latencies in raster plot)
% vecRestrictRange = [0 inf];%do we want to restrict the peak detection to for example the time during stimulus? Then put [0 1] here.
% boolDirectQuantile = false;%if true; uses the empirical null distribution rather than the Gumbel approximation. Note that in this case the accuracy of your p-value is limited by the # of resamplings
% dblJitterSize = 2; %scalar value, sets the temporal jitter window relative to dblUseMaxDur (default: 2). Note that a value of "2" therefore means that the last data point that can be used is at t=last_onset + dblUseMaxDur*3
% % boolStitch = false; %boolean, stitches data in the case of heterogeneous inter-event durations (default: true)
% 
% %then run ZETA with those parameters
% hTic2=tic;
% dblElapsedTime2=toc(hTic2);
% % fprintf("\nSpecified parameters (elapsed time: %.2f s):\nzeta-test p-value: %f\nt-test p-value: %f\n",dblElapsedTime2,dblZetaP,sZETA.dblMeanP)
% 
% 
% %then run ZETA with those parameters
% % hTic2=tic;
% % [dblZetaP,sZETA,sRate,sLatencies] = zetatest(vecSpikeTimes1,matEventTimes,dblUseMaxDur,intPlot,vecRestrictRange,boolDirectQuantile);
% % dblElapsedTime2=toc(hTic2);
% % fprintf("\nSpecified parameters (elapsed time: %.2f s):\nzeta-test p-value: %f\nt-test p-value: %f\n",dblElapsedTime2,dblZetaP,sZETA.dblMeanP)
% 

end