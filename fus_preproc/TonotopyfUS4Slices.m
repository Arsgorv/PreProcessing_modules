%% Plot tonotopy 4 slices altogether:

%% Loadind all slices

Slices = {'A', 'B', 'C', 'D'};
datapath = sessions{1};
results = struct();

for iSlice = 1:length(Slices)
    fus_file = dir(strcat(datapath, '/fUS/RP_data_*slice_', Slices{iSlice}, '.mat'));
    load(fullfile(datapath, 'fUS', fus_file.name))
    load(fullfile(datapath, 'Master_sync.mat'))
    results.(Slices{iSlice}).cat_tsd = cat_tsd;
end
%% Computing preferred response for each slice

WindowStart = 2;% n seconds post stim onset for PSTH window
WindowEnd = 4; % n seconds post stim onset for PSTH window
TypeBaseline = '10_trials'; % Baseline computed only based on Current_trial or on 10_trials
FREQS =[200 400 800 1600 3200 6400 12800 ];% Frequencis we want to display
FREQSAll =[200 400 800 1600 3200 6400 12800 25600]; % all frequencies all tonotopies (to have same scale across recordings)


for iSlice = 1:length(Slices)
    cat_tsd = results.(Slices{iSlice}).cat_tsd;
    [idxMap, respAmp] = TonotopyPreferredFreq(cat_tsd, Epochs, WindowStart, WindowEnd, FREQS, FREQSAll, TypeBaseline);
    results.(Slices{iSlice}).idxMap  = idxMap;
    results.(Slices{iSlice}).respAmp = respAmp;
    
end 

%% Thresholding preferred responses map based on maximal responses in the auditory cortex

SliceForThreshold = 'B';
thresholdRatio = 0.1;   % 5 

cat_tsd = results.(SliceForThreshold).cat_tsd;
fUSData =  Data(cat_tsd.data);
fUSDataReshaped = reshape(fUSData', cat_tsd.Nx, cat_tsd.Ny, size(fUSData, 1)');

figure(12);
imagesc(mean(fUSDataReshaped,3))
axis image
colormap hot
title('Draw mask: ACx ROI to fix threshold');

ACmask = roipoly;
respAmp = results.B.respAmp;
respMaxPix = max(respAmp, [], 3);
respAC = respMaxPix(ACmask);
RefResponse = prctile(respAC(:), 95);

Threshold = thresholdRatio * RefResponse;

for iSlice = 1:4
    
    respAmp = results.(Slices{iSlice}).respAmp;
    idxMap  = results.(Slices{iSlice}).idxMap;
    
    respMaxPix = max(respAmp, [], 3);
    
    idxMapThresh = idxMap;
    idxMapThresh(respMaxPix < Threshold) = 0;
    
    results.(Slices{iSlice}).idxMapThresh = idxMapThresh;
end

%% Save Best Frequency map and mean fUS image in a struct 

MatBestFreq = struct();
for iSlice=1:4
    MatBestFreq.(Slices{iSlice}).BestFrequencyMap =results.(Slices{iSlice}).idxMap;
    MatBestFreq.(Slices{iSlice}).BestFrequencyMapThesholded =results.(Slices{iSlice}).idxMapThresh;
    
    cat_tsd = results.(Slices{iSlice}).cat_tsd;
    
    fUSData =  Data(cat_tsd.data);
    fUSDataReshaped = reshape(fUSData', cat_tsd.Nx, cat_tsd.Ny, size(fUSData, 1)');
    
    MatBestFreq.(Slices{iSlice}).MeanImage = mean(fUSDataReshaped,3);
end


save(fullfile(datapath, 'MatBestFreq.mat'),'MatBestFreq')

%% Visualization 
figure('units','normalized','outerposition',[0.05 0.2 0.9 0.6]);
parts = strsplit(datapath, filesep); 
Ferret = parts{end-1}; 
sessionName = parts{end}; 
sgtitle(strcat(Ferret, ' ', sessionName));
cmapTonotopy = hsv(length(FREQSAll+1));
cmapTonotopy(1,:) = [0 0 0];


for iSlice=1:length(Slices)
    ax1 = subplot(2, length(Slices),iSlice);
    cat_tsd = results.(Slices{iSlice}).cat_tsd;
    
    fUSData =  Data(cat_tsd.data);
    fUSDataReshaped = reshape(fUSData', cat_tsd.Nx, cat_tsd.Ny, size(fUSData, 1)');
    
    imagesc(mean(fUSDataReshaped,3));

   % axis image ij;
    colormap(ax1, hot);
    colorbar;
%     roiFields = { ...
%         'MEG', ...  
%         'AEG', ...
%         'PEG', ...
%         'VP', ...
%         'Hpc', ...        % Hippocampus (legacy: Hpc)
%         'Thal', ...       % Thalamus
%         'tissue', ...     % Tissue (legacy: tissue)
%         'vessel', ...     % Vessel (legacy: vessel)
%         'far_out'};   
%     hold on;
%     if exist('masks', 'var')
%         for r = 1:numel(roiFields)
%             fName = roiFields{r};
%             if isfield(masks, fName) && ~isempty(masks.(fName))
%                 thisMask = masks.(fName);
%                 if isequal(size(thisMask), size(fUSDataReshaped(:,:,1)))
%                     contour(thisMask, [0.5 0.5], 'Color', 'r', 'LineWidth', 0.5);
%                 end
%             end
%         end
%     end
    title(strcat('Slice', Slices{iSlice}));
    ax2 = subplot(2, length(Slices),iSlice + length(Slices));
    idxMapThresh = results.(Slices{iSlice}).idxMapThresh;
    imagesc(idxMapThresh);
    %axis image ij;
    nFAll = numel(FREQSAll);
    colormap(ax2, cmapTonotopy);
    caxis([0 nFAll]);
    colorbar('Ticks', 0:nFAll, 'TickLabels', string([0, FREQSAll]))
    linkaxes([ax1 ax2],'xy');

end 

saveas(gcf, fullfile(datapath, 'figures/TonotopyPreferredResponse.png'))

 

function [idxMap, respAmp] = TonotopyPreferredFreq(cat_tsd, Epochs, WindowStart, WindowEnd, FREQS, FREQSAll, TypeBaseline)
    
    if size(FREQS) == size(FREQSAll)
        if ~isfield(Epochs, sprintf('stim_%d', FREQS(end)))
            FREQS(end) = [];
        end
    end
    nPreDurationSeconds = 2;
    nPostDurationSeconds = 5.5;
    postStimDur =nPostDurationSeconds * 10000; 
    preStimDur =nPreDurationSeconds * 10000; 

    nF = numel(FREQS);
    PSTH = cell(nF,1);
    timePSTH = cell(nF,1);

    for iF = 1:nF
        
        % Extract all stim corresponding to each frequency
        
        stim_epochs = Epochs.(sprintf('stim_%d', FREQS(iF)));

        tStart = Start(stim_epochs);
        tEnd   = End(stim_epochs);
        
        nTrials = length(tStart); 
        for iTrial = 1:nTrials
            
            trial_epoch = intervalSet( tStart(iTrial) - preStimDur,tEnd(iTrial) + postStimDur);
            restr = Restrict(cat_tsd.data, trial_epoch);
            trial_data = Data(restr);
            trial_time = Range(restr);
            
            trial_data = reshape(trial_data', cat_tsd.Nx, cat_tsd.Ny, []);
            
            TrialsData{iTrial} = trial_data;
            TrialsTime{iTrial} = trial_time - trial_time(1);
        end

        %Baseline Correction for each pixel and each trial:
        % either based on
        % Current_trial pre stim silence
        % 10_trials: baseline computed on pre-stim silences averages across the 5
            %precedings and following trials.
        
        AllStimStarts = Start(Epochs.stim_all);
        TrialsBaselineCorr = {};
        for iTrial = 1:length(TrialsData)
            trial = TrialsData{iTrial};
            if strcmp(TypeBaseline, 'Current_trial')
                baseline = mean(trial(:,:,1:2*nPreDurationSeconds), 3);
                TrialsBaselineCorr{iTrial} = (trial - baseline)./baseline;
            elseif strcmp(TypeBaseline, '10_trials')
                stim_begins_trial = tStart(iTrial);
                idxTrialWithinSession =  find(AllStimStarts == stim_begins_trial);

                idxStartBaseline = max(1, idxTrialWithinSession - 5);
                idxEndBaseline   = min(length(AllStimStarts), idxTrialWithinSession + 5);

                baselineSegments = [];
                AllTrialStarts = Start(Epochs.trial_all);
                for jTrial = idxStartBaseline:idxEndBaseline
                    trialStart = AllTrialStarts(jTrial);

                    SilencePreTrialEpoch = intervalSet(trialStart, trialStart+preStimDur);

                    restr = Restrict(cat_tsd.data, SilencePreTrialEpoch);
                    data  = Data(restr);
                    data = reshape(data', cat_tsd.Nx, cat_tsd.Ny, []);
                    baselineSegments = cat(3, baselineSegments, data);
                end
                baseline = mean(baselineSegments, 3);
                TrialsBaselineCorr{iTrial} = (trial - baseline) ./ baseline;
            end
                
        end
        
        % Compute PSTH for each Freq
        
         %checkin that all trials have same length
        segLengths = cellfun(@(x) size(x,3), TrialsBaselineCorr);
        targetLen = mode(segLengths);
        nSeg = numel(TrialsBaselineCorr);
        stack = nan(cat_tsd.Nx, cat_tsd.Ny, targetLen, nSeg);

        for i = 1:nSeg
            seg = TrialsBaselineCorr{i};
            L = size(seg,3);

            if L >= targetLen
                stack(:,:,:,i) = seg(:,:,1:targetLen);
            end
        end
        PSTH{iF} = mean(stack, 4, 'omitnan'); % Nx × Ny × time
        timePSTH{iF} = TrialsTime{1} - TrialsTime{1}(1);
    end
    
    
    % Compute Preferred frequency index based on PSTH and response amplitude
    Nx = cat_tsd.Nx;
    Ny = cat_tsd.Ny;
    tonoMap = nan(Nx, Ny);
    respAmp = nan(Nx, Ny, nF); % amplitude max (optionnel)
    for ix = 1:Nx
        for iy = 1:Ny
            bestResp = -Inf;
            bestFreq = NaN;
            for iF = 1:nF
                t = timePSTH{iF}/10000;
                psth = squeeze(PSTH{iF}(ix, iy, :));

                stimStarts = nPreDurationSeconds;
                idxWin = t >= stimStarts + WindowStart & t <= (stimStarts + WindowEnd);
                respVal = mean(psth(idxWin), 'omitnan');

                if respVal > bestResp
                    bestResp = respVal;
                    bestFreq = FREQS(iF);
                end
                respAmp(ix, iy, iF) = respVal;
            end
            tonoMap(ix, iy) = bestFreq;
        end
    end
    % discrete index ranging from 1:AllnF
    [~, idxMap] = ismember(tonoMap, FREQSAll);  
end

