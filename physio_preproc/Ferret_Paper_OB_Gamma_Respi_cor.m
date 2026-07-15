
clear all
%% make data
MakeRespi_ForSession_Ferret


%% beautiful example
cd('/media/nas8/OB_ferret_AG_BM/Shropshire/head-fixed/20241211_TORCs')
load('SleepScoring_OBGamma.mat', 'RespRate_tsd', 'RespRateVar_tsd', 'CleanStates' ,'SmoothGamma')

X = Data(RespRate_tsd);
Y = log10(Data(Restrict(SmoothGamma , RespRate_tsd)));
Fs = 1250;              % sampling rate
binSize = 10;           % seconds
binSamples = Fs * binSize;

nBins = floor(length(X) / binSamples);

X_binned = zeros(nBins,1);
Y_binned = zeros(nBins,1);
for i = 1:nBins
    idx = (i-1)*binSamples + (1:binSamples);
    X_binned(i) = mean(X(idx));   % average per 10s bin
    Y_binned(i) = mean(Y(idx));   % average per 10s bin
end


figure
PlotCorrelations_BM(X_binned , Y_binned , 'marker_size' , 5)
xlabel('Respiratory rate (Hz)'), ylabel('OB gamma power (log)')
axis square, xlim([.4 1.35]), ylim([2.05 2.6])
makepretty_BM2





%% cross-corr


%%
cd('/media/nas7/React_Passive_AG/OBG/Labneh/head-fixed/20230227/')
load('SleepScoring_OBGamma.mat', 'CleanStates', 'SmoothGamma' , 'RespRate_tsd')

load('LFPData/LFP26.mat')
FilGamma = FilterLFP(LFP,[40 60],1024);
tEnveloppeGamma = tsd(Range(LFP), abs(hilbert(Data(FilGamma))) ); 
smootime=.06;
SmoothGamma = tsd(Range(tEnveloppeGamma), runmean(Data(tEnveloppeGamma), ...
    ceil(smootime/median(diff(Range(tEnveloppeGamma,'s'))))));


FilULow = FilterLFP(LFP,[.1 1],1024); 
tEnveloppeULow = tsd(Range(LFP), abs(hilbert(Data(FilULow))) ); 
smootime=.006;
SmoothULow = tsd(Range(tEnveloppeULow), runmean(Data(tEnveloppeULow), ...
    ceil(smootime/median(diff(Range(tEnveloppeULow,'s'))))));


[c_all,lags] = xcorr(Data(SmoothULow) , Data(Restrict(SmoothGamma,SmoothULow)) , 3750);
SmoothULow_Wake = Restrict(SmoothULow , Wake);
[c_wake,lags] = xcorr(Data(SmoothULow_Wake) , Data(Restrict(SmoothGamma,SmoothULow_Wake)) , 3750 , 'biased');
SmoothULow_NREM = Restrict(SmoothULow , SWSEpoch);
[c_nrem,lags] = xcorr(Data(SmoothULow_NREM) , Data(Restrict(SmoothGamma,SmoothULow_NREM)) , 3750 , 'biased');
SmoothULow_REM = Restrict(SmoothULow , REMEpoch);
[c_rem,lags] = xcorr(Data(SmoothULow_REM) , Data(Restrict(SmoothGamma,SmoothULow_REM)) , 3750 , 'biased');


figure
subplot(121)
plot(linspace(-3,3,7501) , c_wake , 'b' , 'LineWidth' , 2)
hold on
plot(linspace(-3,3,7501) , c_nrem , 'r' , 'LineWidth' , 2)
plot(linspace(-3,3,7501) , c_rem , 'g' , 'LineWidth' , 2)
vline(0,'--r')
xlabel('lag (s)'), ylabel('Corr values (a.u.)'), xlim([-3 3])
box off

SmoothULow_Wake = Restrict(SmoothULow , Wake);
[c_wake,lags] = xcorr(zscore(Data(SmoothULow_Wake)) , zscore(Data(Restrict(SmoothGamma,SmoothULow_Wake))) , 3750 , 'biased');
SmoothULow_NREM = Restrict(SmoothULow , SWSEpoch);
[c_nrem,lags] = xcorr(zscore(Data(SmoothULow_NREM)) , zscore(Data(Restrict(SmoothGamma,SmoothULow_NREM))) , 3750 , 'biased');
SmoothULow_REM = Restrict(SmoothULow , REMEpoch);
[c_rem,lags] = xcorr(zscore(Data(SmoothULow_REM)) , zscore(Data(Restrict(SmoothGamma,SmoothULow_REM))) , 3750 , 'biased');

subplot(122)
plot(linspace(-3,3,7501) , c_wake , 'b' , 'LineWidth' , 2)
hold on
plot(linspace(-3,3,7501) , c_nrem , 'r' , 'LineWidth' , 2)
plot(linspace(-3,3,7501) , c_rem , 'g' , 'LineWidth' , 2)
vline(0,'--r')
xlabel('lag (s)'), ylabel('Corr values (a.u.)'), xlim([-3 3])
box off



%% tools
figure
plot(Range(SmoothULow,'s') , Data(SmoothULow))
hold on
plot(Range(SmoothGamma,'s') , Data(SmoothGamma))
xlim([8249 8259])
