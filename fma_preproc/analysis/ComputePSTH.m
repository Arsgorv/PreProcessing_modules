function [meanPSTH, semPSTH, trialsPSTH, meanZ, semZ, trialsPSTH_z] = ComputePSTH(spikeTimes, onsetTimes, timeEdges)
% Helper function to compute PSTH over a set of onsetTimes
% Input:
% spikeTimes: array of spike times for a single neuron (in seconds)
% onsetTimes: array of stimulus onset times
% timeEdges: edges for histogram bins
% Output:
% meanPSTH   : (1 x nBins)
% semPSTH    : (1 x nBins)
% trialsPSTH : (nTrials x nBins), trial-by-trial rates

if nargin < 3
    error('ComputePSTH requires spikeTimes, onsetTimes, and edges.');
end

nBins   = length(timeEdges)-1;
nOnsets = length(onsetTimes);
trialsPSTH = zeros(nOnsets,nBins);
binWidth   = diff(timeEdges(1:2)); % assuming uniform bin size

for iTr=1:nOnsets
    relTimes = spikeTimes - onsetTimes(iTr);
    inWin    = (relTimes>=timeEdges(1) & relTimes<=timeEdges(end));
    theseSpk = relTimes(inWin);
    counts   = histcounts(theseSpk, timeEdges);
    trialsPSTH(iTr,:) = counts / binWidth; % spikes/s
end

meanPSTH = mean(trialsPSTH,1);
stdPSTH  = std(trialsPSTH,0,1);
semPSTH  = stdPSTH / sqrt(size(trialsPSTH, 1));

% Z-score
idxBL = (timeEdges>=-1 & timeEdges<0);
baseline = trialsPSTH(:,idxBL);
meanBL = mean(baseline,'all'); %mean across trials and across -1-0 window
stdBL = std(baseline,0,'all');
trialsPSTH_z = (trialsPSTH - meanBL) ./ stdBL;

meanZ = mean(trialsPSTH_z,1);
stdZ  = std(trialsPSTH_z,0,1);
semZ = stdZ/sqrt(nOnsets);



end
