function [Times,Template,HeartRate,GoodEpoch]=DetectHeartBeats_EmbReact_BM(LFP,NoiseEpoch,Options,plo)
%
% This function returns the time of individual heartbeats
%
% INPUT
%   LFP : LFP of EKG recording
%
% OUTPUT
%   Times : times of beats
%
%
% Author: Baptiste Maheo, changed compared to DetectHeartBeats_EmbReact_SB on line 97
%

%% Suggested input options
try
    Options.TemplateThreshStd;
catch
    Options.TemplateThreshStd=3;
end
try
    Options.BeatThreshStd;
catch
    Options.BeatThreshStd=1.5;
end


%% Important parameters
FilterBand = [30 200];
TemplateDuration = 10; %ms
SnipDur = 5; %s
MinBeatSep = 1/Options.MaxFreq; % ms, so heart rate cannot be above 20Hz
TemplateThreshStd = Options.TemplateThreshStd; % value used to threshold to get the template for convolution
BeatThreshStd = Options.BeatThreshStd; % value used to threshold to final beats


%% init
TotalEpoch = intervalSet(0,max(Range(LFP)));
LFPtemp = Restrict(LFP,TotalEpoch-NoiseEpoch);
Filt_LFP_High = FilterLFP(LFPtemp,FilterBand,1024);
litEpoch = intervalSet(min(Range(Filt_LFP_High)), min(min(Range(Filt_LFP_High))+300*1e4,max(Range(Filt_LFP_High))));
x_data = Range(Restrict(Filt_LFP_High,litEpoch));
y_data = Data(Restrict(Filt_LFP_High,litEpoch));


%% convolve and find match
xmx = findpeaks(y_data, std(y_data) * TemplateThreshStd);
[MR,~] = PlotRipRaw(Filt_LFP_High, x_data(xmx.loc)/1e4, TemplateDuration,1,1);
close

Template = MR(:,2);
ConvSigForNoise = tsd(Range(Filt_LFP_High), conv(Data(Filt_LFP_High),Template,'same') .* conv(Data(Filt_LFP_High), Template, 'same'));
STDLim=std(Data(ConvSigForNoise));
ConvSig = tsd(Range(Filt_LFP_High), conv(Data(Filt_LFP_High),Template,'same') .* conv(Data(Filt_LFP_High), Template, 'same'));

NumSnips = ceil(max(Range(Filt_LFP_High,'s'))/SnipDur);
Times=[]; Vals=[];

for i=1:NumSnips
    LocEpoch = intervalSet((i-1)*SnipDur*1e4,(i)*SnipDur*1e4);
    LocalSigConv = Restrict(ConvSig,LocEpoch);
    LocalSig = Restrict(Filt_LFP_High,LocEpoch);
    if not(isempty(Data(LocalSigConv)))
        x_data = Range(LocalSigConv);
        y_data = Data(LocalSig);
        xmx = findpeaks(Data(LocalSigConv),STDLim*BeatThreshStd);
%         xmx.loc = xmx.loc(y_data(xmx.loc)>std(y_data)*BeatThreshStd);
        Times = [Times;x_data(xmx.loc)];
        Vals = [Vals;y_data(xmx.loc)];
    end
end


%% Get rid of too close peaks
TimeDiff = diff(Times);
problempts = find(TimeDiff<MinBeatSep*1e4);
while (not(isempty(problempts)))
    ToExclude = [];
    for k=1:length(problempts)
        ToCompare = [Vals(problempts(k)),Vals(problempts(k)+1)];
        [~, ind] = min(ToCompare);
        ToExclude(k) = problempts(k)+ind-1;
    end
    Times(ToExclude) = [];
    Vals(ToExclude) = [];
    TimeDiff = diff(Times);
    problempts = find(TimeDiff<MinBeatSep*1e4);
end



%% BM change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% cleaning
data = 1./diff(Times/1e4);
% Parameters
baseWindow = 50;
maxWindow = 200;
outlierThreshold = 0.2;  % Max fraction of outliers per window
zscoreThresh = 2;        % Z-score threshold to detect outliers

n = length(data);
cleanData = data;
outlierMask = false(1, n);

% Sliding window
for i = 1:n
    for w = baseWindow:10:maxWindow
        left = max(1, i - floor(w / 2));
        right = min(n, i + floor(w / 2));
        idx = left:right;
        window = cleanData(idx);
       
        % Exclude NaNs
        valid = ~isnan(window);
        if sum(valid) < w * 0.5
            continue
        end
       
        winData = window(valid);
        mu = mean(winData);
        sigma = std(winData);
       
        % Z-score of central point
        if sigma > 0
            z = abs((cleanData(i) - mu) / sigma);
            % Compute how many points would be marked as outliers in window
            zAll = abs((winData - mu) / sigma);
            if z > zscoreThresh && mean(zAll > zscoreThresh) < outlierThreshold
                outlierMask(i) = true;
            end
        end
        break  % exit after first valid window
    end
end

% Replace outliers with NaN
cleanData(outlierMask) = NaN;

% Interpolate missing values linearly
cleanData = fillmissing(cleanData, 'linear');


% tsd
HeartRate = tsd(Times(1:end-1),cleanData);
BadEpoch = thresholdIntervals(HeartRate,Options.MinFreq,'Direction','Below'); 
BadEpoch = mergeCloseIntervals(BadEpoch,1*1e4);
GoodEpoch = TotalEpoch-BadEpoch;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLOT 
% modif EC to link two subplots on 20/03/2025 and make starts more visible
if plo
    figure;

    ax1 = subplot(211); % First subplot (LFP + detected heartbeats)
    plot(Range(LFP,'s'), Data(LFP), 'k', 'LineWidth', 1.5), hold on % Black for raw LFP
    plot(Range(Filt_LFP_High,'s'), Data(Filt_LFP_High), 'b', 'LineWidth', 1.5), hold on % Blue for filtered LFP
    plot(Times(1:end)/1e4, Vals(1:end), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r') % Red for detected beats
    xlabel('Time (s)'), ylabel('EKG')
    legend({'Raw LFP', 'Filtered LFP', 'Detected Beats'}, 'Location', 'best')
    grid on

    ax2 = subplot(212); % Second subplot (Heart Rate)
    plot(Range(HeartRate,'s'), Data(HeartRate), 'b', 'LineWidth', 1.5) 
    xlabel('Time (s)'), ylabel('Heart Rate (Hz)')
    legend({'Heart Rate'}, 'Location', 'best')
    grid on

    % Link the x-axes of both plots
    linkaxes([ax1, ax2], 'x'); 
end




end