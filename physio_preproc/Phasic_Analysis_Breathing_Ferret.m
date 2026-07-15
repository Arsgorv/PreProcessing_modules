
% cd('/media/nas8/OB_ferret_AG_BM/Shropshire/head-fixed/20241211_TORCs/')


%% OB gamma on breathing
clear all

cd('/media/nas7/React_Passive_AG/OBG/Labneh/head-fixed/20230227')

load('SleepScoring_OBGamma.mat')
Wake = Wake-TotalNoiseEpoch;

load('LFPData/LFP35.mat')
D_zsc = zscore_sliding(Data(LFP) , 1250*5);
D_zsc_tsd = tsd(Range(LFP) , D_zsc);
LFP_breathing = FilterLFP(D_zsc_tsd , [.3 2]);

load('B_Middle_Spectrum.mat')
B_Sp_tsd = tsd(Spectro{2}*1e4 , Spectro{1});


[P,f,VBinnedPhase,L] = PrefPhaseSpectrum_BM(Restrict(LFP_breathing , Wake) , log10(Data(B_Sp_tsd)) , Range(B_Sp_tsd , 's'),...
Spectro{3} , [.3 2] , 30);
f=figure(1); subplot(3,1,1:2), ylim([15 100]), caxis([2.5 3.9]), colormap jet
[P,f,VBinnedPhase,L] = PrefPhaseSpectrum_BM(Restrict(LFP_breathing , Sleep) , log10(Data(B_Sp_tsd)) , Range(B_Sp_tsd , 's'),...
Spectro{3} , [.3 2] , 30);
f=figure(2); subplot(3,1,1:2), ylim([15 100]), caxis([2.5 3.9]), colormap jet




load('B_Low_Spectrum_26.mat')
B_Sp_tsd = tsd(Spectro{2}*1e4 , Spectro{1});


[P,f,VBinnedPhase,L] = PrefPhaseSpectrum_BM(Restrict(LFP_breathing , Wake) , log10(Data(B_Sp_tsd)) , Range(B_Sp_tsd , 's'),...
Spectro{3} , [.3 2] , 30);
f=figure(3); subplot(3,1,1:2), ylim([0 20]), caxis([4 5]), colormap jet
[P,f,VBinnedPhase,L] = PrefPhaseSpectrum_BM(Restrict(LFP_breathing , Sleep) , log10(Data(B_Sp_tsd)) , Range(B_Sp_tsd , 's'),...
Spectro{3} , [.3 2] , 30);
f=figure(4); subplot(3,1,1:2), ylim([0 20]), caxis([4 5]), colormap jet




% --- Inputs ---
% LFP_bulb   : tsd, sampled at Fs (e.g. 30 kHz)
% spikes_OB  : ts object with spike times
clear all
cd('/media/nas8/OB_ferret_AG_BM/Shropshire/head-fixed/20241205_TORCs/')

Fs = 30000;             % sampling rate (Hz)
lowCut  = 0.5;          % Hz
highCut = 4;            % Hz

% load('ChannelsToAnalyse/Bulb_deep.mat')
% load(['LFPData/LFP' num2str(channel) '.mat'])

load('LFPData/LFP105.mat')
D_zsc = zscore_sliding(Data(LFP) , 1250*5);
D_zsc_tsd = tsd(Range(LFP) , D_zsc);
LFP_breathing = FilterLFP(D_zsc_tsd , [.3 2]);

% --- 2. Compute instantaneous phase via Hilbert transform ---
analyticSig = hilbert(Data(LFP_breathing));
phaseLFP = angle(analyticSig);        % radians [-pi, pi]

load('SleepScoring_OBGamma.mat', 'spikes_OB','Wake','Sleep' , 'TotalNoiseEpoch')
Wake = Wake-TotalNoiseEpoch;

PhaseTsd = tsd(Range(LFP), phaseLFP);

% --- 3. Get phase values at spike times ---
spk_times_Wake = Range(Restrict(spikes_OB , Wake));         % spike times (1e-4 s units)
spk_times_Sleep = Range(Restrict(spikes_OB , Sleep));         % spike times (1e-4 s units)

spike_phase_wake = Data(Restrict(PhaseTsd , spk_times_Wake));
spike_phase_sleep = Data(Restrict(PhaseTsd , spk_times_Sleep));

% --- 4. Plot polar histogram of spike phases ---
figure
subplot(121)
polarhistogram(spike_phase_wake, 20 , 'FaceColor', [0 0 1], 'FaceAlpha', 0.7)       % 20 bins
title('Wake')

subplot(122)
polarhistogram(spike_phase_sleep, 20 , 'FaceColor', [.3 .3 .3], 'FaceAlpha', 0.7)       % 20 bins
title('Sleep')





