



%% Respi
try
    load('LFPData/LFP105.mat')
catch
    load('LFPData/LFP35.mat')
end
D_zsc = zscore_sliding(Data(LFP) , 1250*5);
D_zsc_tsd = tsd(Range(LFP) , D_zsc);
LFP_breathing = FilterLFP(D_zsc_tsd , [.3 2.5]);
x_f = Data(LFP_breathing);
fs     = 1250;       % sampling rate
t_1e4 = Range(LFP);                 % times in 1e-4 s

% --- Instantaneous phase and breathing rate (Hz) via Hilbert transform ---
hx    = hilbert(x_f);
phi   = unwrap(angle(hx));          % radians
% Instantaneous frequency f = (fs/2pi) * derivative(phase)
dphi  = gradient(phi);              % rad/sample
f_inst = (fs/(2*pi)) * dphi;        % Hz, same length as x

% --- Clean up: valid physiological range & smoothing ---
% Keep only plausible breathing rates (rats ~0.5–15 Hz). Out-of-range -> NaN.
f_inst( (f_inst < 0.3) | (f_inst > 2) ) = NaN;

% Light smoothing (moving median then mean), without distorting edges much
win_med = max(1, round(2 * fs));  % ~200 ms
win_mean = max(1, round(2 * fs)); % ~100 ms
f_med  = movmedian(f_inst, win_med, 'omitnan');
f_smooth = movmean(f_med, win_mean, 'omitnan');

% --- Build Breathing_tsd (timestamps in 1e-4 s; column vectors) ---
RespRate_tsd = tsd(t_1e4, f_smooth);

figure, plot(Range(RespRate_tsd , 's') , Data(RespRate_tsd))

% --- Save to StateEpochSB.mat (append) ---
save('SleepScoring_OBGamma.mat', 'RespRate_tsd', '-append');


%% Respiratory Variability based on breathing rate (RespRateVar_tsd)
% Respiratory variability (std of respiratory rate) in 2 s bins
% Input:  Breathing_tsd (rate in Hz, timestamps in 1e-4 s)
% Output: RespRateVar_tsd (std(Hz) per 2 s bin, timestamps = bin centers in 1e-4 s)

load('SleepScoring_OBGamma.mat', 'RespRate_tsd');

t_1e4 = Range(RespRate_tsd);      % 1e-4 s
rr_hz = Data(RespRate_tsd);       % Hz
t_s   = t_1e4 * 1e-4;              % seconds

bin_dur = 20;                        % seconds
min_samples = 3;                    % require at least this many samples per bin

t_start = 0;                        % start binning at 0 s (session start)
t_end   = max(t_s);                 % end at last sample time
edges   = (t_start:bin_dur:(t_end + bin_dur))';   % column vector of edges
n_bins  = numel(edges) - 1;

RespRateVar = nan(n_bins, 1);
BinCenters  = edges(1:end-1) + bin_dur/2;

% assign each sample to a bin
[~,~,bin_idx] = histcounts(t_s, edges);

for b = 1:n_bins
    idx = (bin_idx == b);
    if any(idx)
        x = rr_hz(idx);
        x = x(isfinite(x));        % keep NaNs as missing
        if numel(x) >= min_samples
            RespRateVar(b) = std(x);   % Hz
        else
            RespRateVar(b) = NaN;
        end
    else
        RespRateVar(b) = NaN;
    end
end

RespRateVar_tsd = tsd((BinCenters * 1e4), RespRateVar);

figure, plot(Range(RespRateVar_tsd , 's') , Data(RespRateVar_tsd))

save('SleepScoring_OBGamma.mat', 'RespRateVar_tsd', '-append');




