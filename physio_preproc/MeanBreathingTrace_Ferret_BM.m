
clear all

cd('/media/nas7/React_Passive_AG/OBG/Labneh/head-fixed/20230227')

load('SleepScoring_OBGamma.mat')
Wake = Wake-TotalNoiseEpoch;

load('LFPData/LFP35.mat')
D_zsc = zscore_sliding(Data(LFP) , 1250*5);
D_zsc_tsd = tsd(Range(LFP) , D_zsc);
Piezo = FilterLFP(D_zsc_tsd , [.3 2]);

% ===========================
% Average Breathing Cycle (Wake vs Sleep)
% ===========================

% --- Parameters ---
fs = 1/median(diff(Range(Piezo,'s')));   % sampling frequency
nPoints = 200;                           % resample each cycle to 200 points
lowpass_cutoff = 20;                     % Hz, low-pass filter for piezo

% --- Step 1: Filter piezo ---
[b,a] = butter(4, lowpass_cutoff/(fs/2), 'low');
PiezoFilt = tsd(Range(Piezo), filtfilt(b,a,Data(Piezo)));

% --- Step 2: Restrict signal to Wake and Sleep ---
PiezoWake  = Restrict(PiezoFilt, Wake);
PiezoSleep = Restrict(PiezoFilt, Sleep);

% --- Step 3: Detect inhalation troughs ---
% Wake
tW   = Range(PiezoWake,'s');
sigW = Data(PiezoWake);
locsW = findpeaks(-sigW);   % indices of troughs
% Filter out shallow troughs
thrW = -std(sigW)/2;
validW = sigW(locsW.loc) < thrW;
locsW.loc = locsW.loc(validW);
InhTimesWake = tW(locsW.loc);

% Sleep
tS   = Range(PiezoSleep,'s');
sigS = Data(PiezoSleep);
locsS = findpeaks(-sigS);
thrS = -std(sigS)/2;
validS = sigS(locsS.loc) < thrS;
locsS.loc = locsS.loc(validS);
InhTimesSleep = tS(locsS.loc);

% --- Step 4: Extract and normalize cycles ---
cyclesWake  = [];
for i = 1:length(InhTimesWake)-1
    idx = tW >= InhTimesWake(i) & tW < InhTimesWake(i+1);
    if sum(idx) > 5
        segt = linspace(0,1,sum(idx));
        segd = sigW(idx);
        segn = interp1(segt, segd, linspace(0,1,nPoints));
        cyclesWake = [cyclesWake; segn];
    end
end

cyclesSleep = [];
for i = 1:length(InhTimesSleep)-1
    idx = tS >= InhTimesSleep(i) & tS < InhTimesSleep(i+1);
    if sum(idx) > 5
        segt = linspace(0,1,sum(idx));
        segd = sigS(idx);
        segn = interp1(segt, segd, linspace(0,1,nPoints));
        cyclesSleep = [cyclesSleep; segn];
    end
end

% --- Step 5: Average waveform ---
Mwake  = mean(cyclesWake,1);
Msleep = mean(cyclesSleep,1);
SEwake  = std(cyclesWake,[],1)/sqrt(size(cyclesWake,1));
SEsleep = std(cyclesSleep,[],1)/sqrt(size(cyclesSleep,1));

% --- Step 6: Plot ---
figure('Color','w'); hold on
x = linspace(0,1,nPoints);

% Shaded error for Wake
fill([x fliplr(x)], [Mwake-SEwake fliplr(Mwake+SEwake)], ...
    [0.5 0.5 1], 'EdgeColor','none','FaceAlpha',0.3);
plot(x, Mwake,'b','LineWidth',2);

% Shaded error for Sleep
fill([x fliplr(x)], [Msleep-SEsleep fliplr(Msleep+SEsleep)], ...
    [1 0.5 0.5], 'EdgeColor','none','FaceAlpha',0.3);
plot(x, Msleep,'r','LineWidth',2);

xlabel('Normalized breathing cycle (0=inhale trough → 1=next trough)');
ylabel('Piezo amplitude (a.u.)');
legend({'Wake','Sleep'});
title('Average breathing cycle waveform');







