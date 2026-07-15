
function MakeHeartRateForSession_BM(varargin)

% default values for mice
for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        error(['Parameter ' num2str(i+2) ' is not a property.']);
    end
    switch lower(varargin{i})
        case 'lowest_freq'
            lowest_freq = varargin{i+1};
        case 'highest_freq'
            highest_freq = varargin{i+1};
        case 'species'
            species = varargin{i+1};            
    end
end

if species == "ferret"
    lowest_freq = 2; highest_freq = 6;
elseif species == "mouse"
    lowest_freq = 5; highest_freq = 14;
    elseif species == "rat"
    lowest_freq = 3; highest_freq = 10;
end

if ~exist('lowest_freq','var')
    lowest_freq = 5;
end
if ~exist('highest_freq','var')
    highest_freq = 14;
end


clear TTLInfo Behav EKG channel
close all
Options.TemplateThreshStd=3;
Options.BeatThreshStd=0.05;
Options.MaxFreq = highest_freq; 
Options.MinFreq = lowest_freq; 

load('ChannelsToAnalyse/EKG.mat')
load(['LFPData/LFP',num2str(channel),'.mat'])
try
    load('SleepScoring_OBGamma.mat', 'TotalNoiseEpoch')
catch
    load('StateEpochSB.mat', 'TotalNoiseEpoch')
end
load('ExpeInfo.mat')
if ExpeInfo.SleepSession==0
    load('behavResources.mat')
    try,  TTLInfo;
        NoiseEpoch=or(TotalNoiseEpoch,intervalSet(Start(TTLInfo.StimEpoch),Start(TTLInfo.StimEpoch)+2*1e4));
    catch
        NoiseEpoch=TotalNoiseEpoch;
    end
else
    NoiseEpoch=TotalNoiseEpoch;
end
[Times,Template,HeartRate,GoodEpoch] = DetectHeartBeats_EmbReact_BM(LFP,NoiseEpoch,Options,1);
EKG.HBTimes=ts(Times);
EKG.HBShape=Template;
EKG.DetectionOptions=Options;
EKG.HBRate=HeartRate;
EKG.GoodEpoch=GoodEpoch;

save('HeartBeatInfo.mat','EKG')
saveas(1,'EKGCheck.fig'),
saveas(1,'EKGCheck.png')
% close all
clear EKG NoiseEpoch TotalNoiseEpoch TTLInfo LFP EKG HearRate Template Times

close
