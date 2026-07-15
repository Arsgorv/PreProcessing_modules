

function Sniff = Detect_Sniff_From_Piezzo(piezzo)

% Input: piezzo  (tsd)  % thoracic cage signal

% Get time and data
t = Range(piezzo);                 % in 1e-4 s units
x = Data(piezzo);

% Sampling rate (Hz)
dt_sec = median(diff(t)) / 1e4;
Fs = 1 / dt_sec;

% zscore with a 10s window
x = zscore_sliding(x , Fs*5);

% Band-pass filter 0.3–2 Hz
FilLFP = FilterLFP(tsd(t , x) , [.3 2] , 1024);
x_f = Data(FilLFP);

% Hilbert for amplitude (envelope) to reject tiny cycles
z   = hilbert(x_f);
env = abs(z);

% Rising zero-crossings (inspiration onsets)
zc = find(x_f(1:end-1) <= 0 & x_f(2:end) > 0) + 1;

% Amplitude gate (reject weak/ noisy crossings)
% env_thr = prctile(env, 30);        % adjust (e.g., 20–40) if needed
% zc = zc(env(zc) >= env_thr);

% Refractory period so detections stay within 0.3–2 Hz range
% min_ibi_sec = 0.30;                % minimum 0.30 s between onsets
% min_ibi_smp = round(min_ibi_sec * Fs);
% keep = zc([true; diff(zc) >= min_ibi_smp]);

% Event times (keep column vector) and ts object
insp_times = t(zc);
insp_times = insp_times(:);
Sniff = ts(insp_times);

% (optional) quick check plot
figure; plot(t/1e4, x_f); hold on
plot(insp_times/1e4, x_f(zc), 'ro'); xlabel('Time (s)'); ylabel('Piezzo (filtered)');
title('Detected inspiration onsets (0.3–2 Hz)')
xlim([1000 1010])








