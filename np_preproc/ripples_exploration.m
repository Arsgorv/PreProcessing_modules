%% ===================================================================
%  ripples_exploration.m  -  ferret hippocampal SWR exploration
%  -------------------------------------------------------------------
%  Workflow (matches the 3 requested stages):
%    A. setup
%    B. localize CA1 layers on the probe        (layers currently unknown)
%    C. study MANUAL candidates + measure band   (stage 1)
%    D. define detection parameters              (stage 2)
%    E. automatic detection                      (stage 3a)
%    F. study AUTOMATIC candidates + compare      (stage 3b)
%
%  Ferret-specific notes / departures from rodent defaults:
%   - Ripple peak frequency scales DOWN with brain size. Ferret sits
%     between rat (~140-200 Hz) and cat/macaque (~110-125 Hz), so we
%     EXPECT a peak around ~100-150 Hz, likely below the rat band.
%     => Do NOT hard-code [120 220]. Measure the endogenous peak (cell C2)
%        and set rippleBand from it.
%   - A true CA1 SWR is two coupled events: the ripple (pyramidal layer)
%     AND a sharp wave (sink in str. radiatum). We require both + a CSD
%     dipole to reject gamma / filtered spikes / artifacts.
%   - State gating (NREM/EMG/accelero) is not available YET. Detection
%     runs on the whole recording; an epoch hook (detectEpochs) is left
%     in place so it can be restricted later without touching the core.
%
%  Ref: Liu et al. 2022 Nat Commun (consensus on SWR detection),
%       doi:10.1038/s41467-022-33536-x ; Buzsaki/Logothetis/Singer 2013
%       Neuron, doi:10.1016/j.neuron.2013.10.002 (frequency scaling).
%
%  NOTE: written for MATLAB R2018b. Avoids xline/yline. Uses your
%  existing idioms (butter+filtfilt, hilbert envelope, memmap, MAD z).
%  Untested against your data - run cell by cell and check the sanity
%  figures flagged with "SANITY".
%% ===================================================================


%% A. Setup -----------------------------------------------------------
basepath = datapath;
datFile  = fullfile(basepath,'continuous.dat');

% Manual overrides (leave [] to auto-read). Use these if the metadata picks
% the wrong stream (e.g. AP 30 kHz instead of the LFP ~2.5 kHz).
fsManual  = [];      % e.g. 2500
nChManual = [];      % e.g. 384
lsbManual = [];      % e.g. 0.195   (uV per bit)

% fs / nChannels / bit-to-uV from CellExplorer session.mat if present, else
% from the OpenEphys structure.oebin (works on sessions without session.mat).
[fs, nCh, lsb, ycoords] = load_session_info(basepath);
if ~isempty(fsManual),  fs  = fsManual;  end
if ~isempty(nChManual), nCh = nChManual; end
if ~isempty(lsbManual), lsb = lsbManual; end

fileInfo      = dir(datFile);
nSamplesTotal = fileInfo.bytes / 2 / nCh;
m = memmapfile(datFile,'Format',{'int16',[nCh nSamplesTotal],'x'});
fprintf('Session: fs=%g Hz, nCh=%d, duration=%.1f s (%.1f min)\n', ...
    fs, nCh, nSamplesTotal/fs, nSamplesTotal/fs/60);

% Channels spanning the hippocampus on the probe (adjust to your probe).
% Exclude bad/dead channels: the SANITY plots showed single-channel CSD
% stripes near ch 45 and 59 - inspect and list them here.
badCh = [45 59];                     % dead/artifact channels
hpcCh = setdiff(1:100, badCh);

% Order hippocampal channels by DEPTH, not by index. Use chanCoords from the
% session if available (returned by load_session_info); else reconstruct from
% the probe geometry (standard contiguous Neuropixels = 2 sites / 20um row).
probeType = 'NP1.0';                        % <-- 'NP1.0' or 'NP2.0' (geometry fallback)
if isempty(ycoords)
    if strcmpi(probeType,'NP2.0'), vpitch = 15; else, vpitch = 20; end
    ycoords = vpitch * floor((0:nCh-1)'/2);
    warning('chanCoords reconstructed from %s geometry. Confirm a smooth ripple-power depth profile.', probeType);
end
[~,depthSort] = sort(ycoords(hpcCh),'ascend');   % stable sort: ties keep index order
chOrder = hpcCh(depthSort);

% MANUAL candidates: list of event peak times (s) you picked by eye.
% You need ~5-15 clear ones to estimate band/duration/amplitude well.
manualEvents = [4730.88 4736.127 4737.554 4741.413 4745.634 4750.015 ...
                1352.200 1396.512 1406.155 1407.278 1411.228 1441.444 1444.770 1453.353];
if any(manualEvents<0 | manualEvents>nSamplesTotal/fs)
    error('manualEvents has times outside [0 %.1f] s for this session - update them.', nSamplesTotal/fs);
end


%% B. Localize CA1 layers (ripple power + sharp-wave CSD) --------------
%  Goal: pick a pyramidal-layer channel (pyrCh, max ripple power) and a
%  str. radiatum channel (radCh, sharp-wave sink) from the depth profile,
%  averaged over the manual events.
%
%  Provisional, deliberately WIDE ferret band for the profile so we don't
%  pre-bias the layer estimate by assuming a rat band.
locBand = [90 180];
locWin  = [-0.3 0.3];
swBand  = [2 40];                 % sharp wave = slow deflection in radiatum

sRel  = round(locWin(1)*fs):round(locWin(2)*fs);
tLoc  = sRel / fs;
evMask = tLoc > -0.03 & tLoc < 0.05;     % around the ripple/SW

[bLoc,aLoc] = butter(3, locBand/(fs/2),'bandpass');
[bSW ,aSW ] = butter(3, swBand /(fs/2),'bandpass');

nC = numel(chOrder);
ripplePowProfile = zeros(nC,1);          % mean peak ripple env per channel
swProfile        = zeros(nC,length(sRel)); % mean low-pass LFP per channel
nUsed = 0;

for e = 1:numel(manualEvents)
    c0  = round(manualEvents(e)*fs);
    idx = c0 + sRel;
    if idx(1) < 1 || idx(end) > nSamplesTotal, continue; end

    xx = double(m.Data.x(chOrder,idx)) * lsb;     % [depth x time]

    % ripple envelope per channel (filtfilt column-wise -> transpose)
    rip   = filtfilt(bLoc,aLoc, xx')';
    envCh = abs(hilbert(rip')');
    ripplePowProfile = ripplePowProfile + max(envCh(:,evMask),[],2);

    % sharp wave (slow) per channel for CSD
    lp = filtfilt(bSW,aSW, xx')';
    swProfile = swProfile + lp;

    nUsed = nUsed + 1;
end
ripplePowProfile = ripplePowProfile / max(1,nUsed);
swProfile        = swProfile / max(1,nUsed);

% CSD of the sharp wave: -d2/dz2 of the low-passed LFP (uniform spacing).
% A current SINK (CA3 input in str. radiatum) shows as a negative LFP
% trough and a CSD sink. VERIFY the sign convention against a known event.
csd = -diff(swProfile,2,1);
csd = [zeros(1,size(csd,2)); csd; zeros(1,size(csd,2))];   % pad to nC

% pyramidal channel = max ripple power
[~,iPyr] = max(ripplePowProfile);
pyrCh = chOrder(iPyr);

% radiatum channel = strongest sharp-wave sink NEAR the pyramidal layer.
% Use the SW LFP (CSD is corrupted by any single bad channel), search a
% window on BOTH sides of pyr (probe orientation unknown) and take the
% largest sink. VERIFY the chosen side against the SANITY depth plot.
% Robust radiatum pick: smooth across channels, take the sharp-wave value
% at the SW time (central window), choose the global SINK (most negative),
% excluding probe ends. The window-around-pyr version kept getting captured
% by the opposite (source) pole, which flips the sharp wave to positive.
swSmooth = movmean(swProfile, 5, 1);
ctrMask  = tLoc > -0.015 & tLoc < 0.015;
swAtPk   = mean(swSmooth(:,ctrMask), 2);
swAtPk(1:5) = NaN; swAtPk(end-4:end) = NaN;     % drop probe ends
[~,iRad] = min(swAtPk);                          % global sharp-wave sink
radCh = chOrder(iRad);

% Manual override: after viewing the SANITY depth plot, set radChManual to
% the channel sitting in the clear sharp-wave SINK, then re-run from cell B.
if exist('radChManual','var') && ~isempty(radChManual)
    radCh = radChManual;
    iRad  = find(chOrder==radCh);
end

fprintf('Localized: pyrCh = %d (depth-rank %d), radCh = %d (depth-rank %d)\n', ...
    pyrCh, iPyr, radCh, iRad);

% SANITY: ripple power should peak in a narrow band of channels (pyr layer);
% the SW/CSD sink should sit a few hundred um below it.
figure('Name','SANITY: layer localization')
subplot(1,3,1)
plot(ripplePowProfile, 1:nC, 'k.-'); axis ij
hold on; plot(ripplePowProfile(iPyr), iPyr, 'ro');
xlabel('peak ripple env'); ylabel('depth rank (1=top)'); title('ripple power')
subplot(1,3,2)
imagesc(tLoc, 1:nC, swProfile); axis xy
xlabel('time (s)'); title('mean SW LFP'); colormap('viridis'); colorbar
subplot(1,3,3)
imagesc(tLoc, 1:nC, csd); axis xy
hold on; plot(0, iRad, 'rx');
xlabel('time (s)'); title('SW CSD (look for sink)'); colorbar


%% C1. Study a MANUAL candidate (single-event multipanel) -------------
%  Pick which manual event to view.
evIdx = 1;
t0    = manualEvents(evIdx);
win   = [-0.5 0.5];

s1 = max(1, floor((t0 + win(1))*fs));
s2 = floor((t0 + win(2))*fs);
idx = s1:s2;
t   = (idx / fs) - t0;

xx  = double(m.Data.x(chOrder,idx)) * lsb;      % depth-ordered
iPyrLoc = find(chOrder==pyrCh);
iRadLoc = find(chOrder==radCh);

% local reference (CAR over hpc channels) - same idea as your draft.
% NOTE: CAR across the ripple-generating layer can attenuate the ripple
% itself; for the cleanest confirmation use CSD (cell shows both).
car  = median(xx,1);
pyr  = xx(iPyrLoc,:) - car;
rad  = xx(iRadLoc,:) - car;

% provisional ripple band for viewing (replaced after cell C2)
if ~exist('rippleBand','var'), rippleBand = locBand; end
[bR,aR] = butter(3, rippleBand/(fs/2),'bandpass');
[bS,aS] = butter(3, swBand    /(fs/2),'bandpass');

ripPyr = filtfilt(bR,aR,pyr);
swRad  = filtfilt(bS,aS,rad);

figure('Name','MANUAL candidate')
subplot(5,1,1); plot(t,pyr,'k'); xlim(win); ylabel('pyr raw uV')
title(sprintf('event %d  t0=%.3f  pyrCh=%d radCh=%d',evIdx,t0,pyrCh,radCh))
subplot(5,1,2); plot(t,ripPyr,'k'); xlim(win); ylabel(sprintf('%d-%d Hz',rippleBand))
subplot(5,1,3); plot(t,swRad,'b'); xlim(win); ylabel('rad SW (2-40Hz)')
subplot(5,1,4)
imagesc(t, 1:nC, filtfilt(bR,aR,(xx-car)')'); axis xy
ylabel('depth'); title('ripple-filtered LFP across depth'); colormap('viridis'); caxis([-10 10])
subplot(5,1,5)
spectrogram(pyr, round(0.08*fs), round(0.075*fs), 1:300, fs, 'yaxis')
ylim([40 300]); caxis([-10 10]); title('pyr spectrogram'); xlabel('time (s)')


%% C2. Measure the endogenous ripple band (event PSD vs baseline) ------
%  Set rippleBand from the data instead of assuming a rodent band.
psdWin   = [-0.3 0.3];
evPSDwin = [-0.04 0.04];          % around the event
basePSD  = [-0.30 -0.12];         % flank baseline

sRel = round(psdWin(1)*fs):round(psdWin(2)*fs);
tt   = sRel/fs;
evM  = tt>=evPSDwin(1) & tt<=evPSDwin(2);
baM  = tt>=basePSD(1)  & tt<=basePSD(2);

nfft = 2^nextpow2(round(0.08*fs));
Pev = []; Pba = [];
for e = 1:numel(manualEvents)
    c0  = round(manualEvents(e)*fs);
    idx = c0 + sRel;
    if idx(1)<1 || idx(end)>nSamplesTotal, continue; end
    xx  = double(m.Data.x(chOrder,idx)) * lsb;
    pyr = xx(find(chOrder==pyrCh),:) - median(xx,1);
    [pe,f] = pwelch(pyr(evM), [], [], nfft, fs);
    [pb,~] = pwelch(pyr(baM), [], [], nfft, fs);
    Pev = [Pev pe]; Pba = [Pba pb];
end
Pev = mean(Pev,2); Pba = mean(Pba,2);

% ratio removes the 1/f background; the ripple shows as a narrow bump.
ratio = Pev ./ (Pba + eps);

% Find the ripple peak as a LOCAL maximum of the ratio, searching ABOVE
% ~100 Hz. Below that the ratio is inflated by sharp-wave / low-gamma power
% and its GLOBAL max sits at the low edge - not the ripple (this was the
% bug that returned a spurious 78 Hz). The true bump here is ~190-200 Hz.
srch  = f>=100 & f<=250;
fAxis = f(srch); fAxis = fAxis(:);
rr    = ratio(srch); rr = rr(:);

% local maxima of rr without relying on a specific findpeaks on the path.
isPeak = false(numel(rr),1);
isPeak(2:end-1) = rr(2:end-1) > rr(1:end-2) & rr(2:end-1) > rr(3:end);
loc = find(isPeak);
if isempty(loc)
    [~,ip] = max(rr);                 % fallback: no interior peak
else
    [~,j] = max(rr(loc)); ip = loc(j);% strongest local peak
end
fPeak = fAxis(ip);

% band around the peak, capped below Nyquist. OVERRIDE if needed.
nyq = fs/2;
rippleBand = round([0.80*fPeak 1.15*fPeak]/5)*5;
rippleBand(2) = min(rippleBand(2), floor(0.9*nyq));
fprintf('Endogenous ripple peak ~ %.0f Hz -> suggested rippleBand = [%d %d] Hz\n', ...
    fPeak, rippleBand(1), rippleBand(2));

figure('Name','SANITY: endogenous band')
subplot(1,2,1)
plot(f,10*log10(Pev),'k'); hold on; plot(f,10*log10(Pba),'Color',[.6 .6 .6])
xlim([0 300]); xlabel('Hz'); ylabel('dB'); legend('event','baseline'); title('PSD')
subplot(1,2,2)
plot(fAxis,rr,'k'); hold on; plot(fPeak,rr(ip),'ro');
xlabel('Hz'); ylabel('event/baseline'); title(sprintf('peak %.0f Hz',fPeak))


%% D. Detection parameters (stage 2) ---------------------------------
%  Set from the manual distributions above. Defaults are starting points;
%  refine using the histograms in cell F1.
detectCh   = pyrCh;               % localized pyramidal layer
swCh       = radCh;               % localized radiatum
refChList  = hpcCh;               % CAR reference set

% rippleBand was set in C2. If you skipped C2, fall back wide:
if ~exist('rippleBand','var'), rippleBand = [100 160]; end

thrHigh    = 5;                   % peak threshold (MAD z)
thrLow     = 2.5;                 % boundary threshold (MAD z)
minDur     = 0.020;              % s  (consensus: >~15-20 ms; check distn)
maxDur     = 0.200;              % s
minInterval= 0.030;              % s  merge events closer than this

% State gating hook: empty = whole recording. Later set to NREM epochs
% as an [nEpoch x 2] matrix of [startSec stopSec].
detectEpochs = [];

fprintf('Detection: ch %d, band [%d %d] Hz, thr %.1f/%.1f, dur %.0f-%.0f ms\n', ...
    detectCh, rippleBand(1), rippleBand(2), thrHigh, thrLow, 1e3*minDur, 1e3*maxDur);


%% E. Automatic detection (stage 3a) ---------------------------------
chunkSec   = 60;
overlapSec = 1;

chunkSamples   = round(chunkSec   * fs);
overlapSamples = round(overlapSec * fs);

[b,a] = butter(3, rippleBand/(fs/2),'bandpass');

channelsToRead = unique([detectCh refChList]);
detectIdx = find(channelsToRead == detectCh);

ripples = [];      % [startSample peakSample stopSample peakZ dur]
chunkStart = 1;

while chunkStart <= nSamplesTotal

    chunkEnd  = min(nSamplesTotal, chunkStart + chunkSamples - 1);
    readStart = max(1, chunkStart - overlapSamples);
    readEnd   = min(nSamplesTotal, chunkEnd + overlapSamples);

    xChunk = double(m.Data.x(channelsToRead,readStart:readEnd)) * lsb;

    x   = xChunk(detectIdx,:);
    ref = median(xChunk(ismember(channelsToRead,refChList),:),1);
    x   = x - ref;

    xf  = filtfilt(b,a,x);
    env = abs(hilbert(xf));

    % robust background. CAVEAT: computed per chunk over ALL states. Once
    % NREM scoring exists, estimate mu/sd over NREM only and apply a fixed
    % threshold, else event RATE becomes state-dependent (Liu et al 2022).
    mu = median(env);
    sd = median(abs(env - mu)) * 1.4826;
    if sd == 0 || isnan(sd)
        chunkStart = chunkEnd + 1; continue
    end
    envZ = (env - mu) / sd;

    aboveHigh = envZ > thrHigh;
    aboveLow  = envZ > thrLow;

    d = diff([0 aboveHigh 0]);
    peakStarts = find(d == 1);
    peakStops  = find(d == -1) - 1;

    keepStart = chunkStart - readStart + 1;
    keepEnd   = chunkEnd   - readStart + 1;

    for i = 1:length(peakStarts)
        p1 = peakStarts(i); p2 = peakStops(i);

        s1 = p1; while s1 > 1 && aboveLow(s1), s1 = s1 - 1; end
        s2 = p2; while s2 < length(aboveLow) && aboveLow(s2), s2 = s2 + 1; end

        [peakVal,iMax] = max(envZ(s1:s2));
        peakSampleLocal = s1 + iMax - 1;

        % assign event to the chunk that owns its peak (avoid double count)
        if peakSampleLocal < keepStart || peakSampleLocal > keepEnd, continue; end

        dur = (s2 - s1 + 1) / fs;
        if dur < minDur || dur > maxDur, continue; end

        ripples = [ripples; ...
            readStart+s1-1, readStart+peakSampleLocal-1, readStart+s2-1, peakVal, dur];
    end

    fprintf('Processed %.1f / %.1f min, %d candidates\n', ...
        chunkEnd/fs/60, nSamplesTotal/fs/60, size(ripples,1));
    chunkStart = chunkEnd + 1;
end

% merge events whose gap < minInterval (consensus: avoid fragmenting one
% ripple into several).
if ~isempty(ripples)
    ripples = sortrows(ripples,1);
    merged  = ripples(1,:);
    for i = 2:size(ripples,1)
        gap = (ripples(i,1) - merged(end,3)) / fs;
        if gap < minInterval
            merged(end,3) = ripples(i,3);                       % extend stop
            if ripples(i,4) > merged(end,4)                     % keep stronger peak
                merged(end,2) = ripples(i,2);
                merged(end,4) = ripples(i,4);
            end
            merged(end,5) = (merged(end,3)-merged(end,1))/fs;   % update dur
        else
            merged = [merged; ripples(i,:)];
        end
    end
    ripples = merged;
end

rippleTimes = ripples(:,2) / fs;

% apply state gating hook (no-op while detectEpochs is empty)
if ~isempty(detectEpochs)
    keep = false(size(rippleTimes));
    for e = 1:size(detectEpochs,1)
        keep = keep | (rippleTimes>=detectEpochs(e,1) & rippleTimes<=detectEpochs(e,2));
    end
    ripples = ripples(keep,:); rippleTimes = rippleTimes(keep);
end

fprintf('\nDetected %d candidate ripples (%.3f /s over whole recording)\n', ...
    size(ripples,1), size(ripples,1)/(nSamplesTotal/fs));


%% F1. Features + sharp-wave/CSD curation (stage 3b) ------------------
%  Per candidate: ripple z, peak frequency, cycles, broadband-contamination
%  ratios, AND a sharp-wave amplitude on radCh. The SW requirement is the
%  main upgrade over the draft - it removes gamma/spike false positives.
win      = [-0.3 0.3];
eventWin = [-0.05 0.05];
baseWin1 = [-0.25 -0.12];
baseWin2 = [ 0.12  0.25];

channelsToRead = unique([detectCh swCh refChList]);
detectIdx = find(channelsToRead == detectCh);
swIdx     = find(channelsToRead == swCh);
refIdx    = find(ismember(channelsToRead,refChList));

[bRipple,aRipple] = butter(3, rippleBand/(fs/2),'bandpass');
[bSWf  ,aSWf]     = butter(3, swBand    /(fs/2),'bandpass');
hiBand = [min(300,floor(0.45*fs)) min(400,floor(0.48*fs))]; % broadband check
[bHigh,aHigh]     = butter(3, hiBand/(fs/2),'bandpass');

sRel = round(win(1)*fs):round(win(2)*fs);
t    = sRel / fs;
eventMask = t>=eventWin(1) & t<=eventWin(2);
baseMask  = (t>=baseWin1(1)&t<=baseWin1(2)) | (t>=baseWin2(1)&t<=baseWin2(2));

N = size(ripples,1);
peakZ=nan(N,1); fPeakEv=nan(N,1); nCycles=nan(N,1);
ratioHigh=nan(N,1); swZ=nan(N,1);

eventSamples = ripples(:,2);
for k = 1:N
    idx = eventSamples(k) + sRel;
    if idx(1)<1 || idx(end)>nSamplesTotal, continue; end

    xx  = double(m.Data.x(channelsToRead,idx)) * lsb;
    car = median(xx(refIdx,:),1);
    pyr = xx(detectIdx,:) - car;
    rad = xx(swIdx,:)     - car;

    fR = filtfilt(bRipple,aRipple,pyr);
    fH = filtfilt(bHigh ,aHigh ,pyr);
    fSWv = filtfilt(bSWf,aSWf,rad);

    env = abs(hilbert(fR));
    mu = median(env(baseMask)); sd = median(abs(env(baseMask)-mu))*1.4826;
    if sd==0 || isnan(sd), continue; end
    peakZ(k) = max((env(eventMask)-mu)/sd);

    % broadband contamination: ripple power should dominate the >300 Hz band
    ratioHigh(k) = mean(fR(eventMask).^2) / (mean(fH(eventMask).^2)+eps);

    % cycles within the event's OWN duration. Counting over a fixed window
    % just measures window_length x frequency (~19 at 195 Hz) and saturates
    % any upper bound - that was rejecting real ripples.
    evStartT = (ripples(k,1)-ripples(k,2))/fs;
    evStopT  = (ripples(k,3)-ripples(k,2))/fs;
    durMask  = t>=evStartT & t<=evStopT;
    y = fR(durMask); nCycles(k) = sum(abs(diff(sign(y)))>0)/2;

    % per-event peak frequency
    [pe,fg] = pwelch(pyr(eventMask),[],[],2^nextpow2(round(0.08*fs)),fs);
    bb = fg>=120 & fg<=240; fg=fg(bb); pe=pe(bb);  % search above the 1/f knee
    [~,ipk]=max(pe); fPeakEv(k)=fg(ipk);           % else argmax sticks at the low edge

    % SHARP WAVE: slow deflection on radCh at the ripple time vs baseline.
    swEnv = abs(hilbert(fSWv));
    muS=median(swEnv(baseMask)); sdS=median(abs(swEnv(baseMask)-muS))*1.4826;
    if sdS>0, swZ(k) = max((swEnv(eventMask)-muS)/sdS); end
end

% curation rule: real ripple = strong narrowband bump + enough cycles +
% not broadband + a co-occurring sharp wave.
good = peakZ>thrHigh & nCycles>=4 & nCycles<=18 & ratioHigh>2 & swZ>2;
fprintf('Candidates: %d  ->  curated (with sharp wave): %d\n', N, sum(good));

ripplesClean    = ripples(good,:);
rippleTimesClean= ripplesClean(:,2)/fs;

% SANITY: feature distributions (report distributions, not means - they are
% log-normal). Look for a unimodal ripple cluster vs an artifact tail.
figure('Name','SANITY: feature distributions')
subplot(2,3,1); histogram(fPeakEv(good)); xlabel('peak freq (Hz)'); title('frequency')
subplot(2,3,2); histogram(1e3*ripples(good,5)); xlabel('duration (ms)'); title('duration')
subplot(2,3,3); histogram(peakZ(good)); xlabel('ripple z'); title('amplitude')
subplot(2,3,4); histogram(nCycles(good)); xlabel('# cycles'); title('cycles')
subplot(2,3,5); histogram(swZ(good)); xlabel('sharp-wave z'); title('sharp wave')
subplot(2,3,6); plot(peakZ,swZ,'.'); hold on; plot(peakZ(good),swZ(good),'.');
xlabel('ripple z'); ylabel('SW z'); title('ripple vs SW')


%% F2. Ripple-triggered average + CSD --------------------------------
%  A clean SWR shows: ripple in pyr, slow trough in radiatum, CSD dipole
%  (sink in radiatum / source near pyr). This is the strongest single check
%  that the curated events are real CA1 SWRs.
plotWin = [-0.15 0.15];
padSec  = 0.25;

sRelFull = round((plotWin(1)-padSec)*fs):round((plotWin(2)+padSec)*fs);
tFull = sRelFull/fs; crop = tFull>=plotWin(1) & tFull<=plotWin(2); t = tFull(crop);

[bR,aR] = butter(3, rippleBand/(fs/2),'bandpass');
[bS,aS] = butter(3, swBand    /(fs/2),'bandpass');

ev = ripplesClean(:,2);
nT = sum(crop);
PyrRaw = zeros(numel(ev),nT); RadSW = PyrRaw; RipEnv = PyrRaw;
sumDepth = zeros(nC,nT); nGood = 0;
for k = 1:numel(ev)
    idx = ev(k)+sRelFull;
    if idx(1)<1 || idx(end)>nSamplesTotal, continue; end
    xx = double(m.Data.x(chOrder,idx)) * lsb; car = median(xx,1);
    pyr = xx(find(chOrder==pyrCh),:)-car;
    rad = xx(find(chOrder==radCh),:)-car;
    radSWfull = filtfilt(bS,aS,rad);
    ripFull   = filtfilt(bR,aR,pyr); envFull = abs(hilbert(ripFull));
    swAll     = filtfilt(bS,aS,(xx-car)')';
    nGood = nGood + 1;
    PyrRaw(nGood,:) = pyr(crop);
    RadSW(nGood,:)  = radSWfull(crop);
    RipEnv(nGood,:) = envFull(crop);
    sumDepth = sumDepth + swAll(:,crop);
end
PyrRaw=PyrRaw(1:nGood,:); RadSW=RadSW(1:nGood,:); RipEnv=RipEnv(1:nGood,:);
meanDepth = sumDepth/nGood;
meanCSD = -diff(meanDepth,2,1); meanCSD = [zeros(1,nT); meanCSD; zeros(1,nT)];

mPyr=mean(PyrRaw,1); sePyr=std(PyrRaw,[],1)/sqrt(nGood);
mSW =mean(RadSW,1);  seSW =std(RadSW,[],1)/sqrt(nGood);
mEnv=mean(RipEnv,1); seEnv=std(RipEnv,[],1)/sqrt(nGood);

figure('Name','F2 ripple-triggered average +/- SEM')
subplot(1,3,1)
fill([t fliplr(t)],[mPyr+sePyr fliplr(mPyr-sePyr)],[.85 .85 .85],'EdgeColor','none','HandleVisibility','off'); hold on
plot(t,mPyr,'k')
fill([t fliplr(t)],[mSW+seSW fliplr(mSW-seSW)],[.8 .8 1],'EdgeColor','none','HandleVisibility','off')
plot(t,mSW,'b')
fill([t fliplr(t)],[mEnv+seEnv fliplr(mEnv-seEnv)],[1 .8 .8],'EdgeColor','none','HandleVisibility','off')
plot(t,mEnv,'r')
xlim(plotWin); legend('pyr raw','rad SW','ripple env')
xlabel('time from peak (s)'); ylabel('\muV'); title(sprintf('mean event \\pm SEM, n=%d',nGood))
subplot(1,3,2); imagesc(t,1:nC,meanDepth); axis xy; xlabel('time (s)'); ylabel('channel (depth)')
title('SW LFP depth'); colormap('viridis'); colorbar
subplot(1,3,3); imagesc(t,1:nC,meanCSD); axis xy; xlabel('time (s)'); ylabel('channel (depth)')
title('SW CSD (blue = sink)'); colorbar


%% F3. Manual-vs-auto comparison (recall + precision spot-check) -------
% Split recall into RAW (pre-curation) vs CURATED so we know whether a miss
% is a detection failure or a curation (sharp-wave gate) failure.
% tolerance is generous: manual marks are imprecise vs the envelope peak.
tol  = 0.10;   % s
rawT = ripples(:,2)/fs;
nRaw = 0; nCln = 0;
fprintf('\nManual event diagnostics (tol %.0f ms):\n', 1e3*tol);
for e = 1:numel(manualEvents)
    [dRaw,iR] = min(abs(rawT - manualEvents(e)));
    [dCln,~]  = min(abs(rippleTimesClean - manualEvents(e)));
    nRaw = nRaw + (dRaw<=tol);
    nCln = nCln + (dCln<=tol);
    % print the curation features of the nearest raw event -> shows WHY it is
    % kept or rejected (which criterion fails).
    fprintf(['  t=%.3f  raw %+5.0f ms | peakZ=%.1f swZ=%.1f ratioHi=%.1f ' ...
             'nCyc=%.0f pass=%d | curated %+6.0f ms\n'], ...
        manualEvents(e), 1e3*(rawT(iR)-manualEvents(e)), ...
        peakZ(iR), swZ(iR), ratioHigh(iR), nCycles(iR), good(iR), 1e3*dCln);
end
fprintf('Recall  raw: %d/%d   curated: %d/%d\n', ...
    nRaw, numel(manualEvents), nCln, numel(manualEvents));

% Precision spot-check: browse a random sample of curated events by eye.
rng(1); nShow = min(20,size(ripplesClean,1));
showIdx = randperm(size(ripplesClean,1), nShow);

plotWin=[-0.15 0.15]; padSec=0.2;
sRelFull=round((plotWin(1)-padSec)*fs):round((plotWin(2)+padSec)*fs);
tFull=sRelFull/fs; crop=tFull>=plotWin(1)&tFull<=plotWin(2); t=tFull(crop);
[bR,aR]=butter(3,rippleBand/(fs/2),'bandpass');
[bS,aS]=butter(3,swBand/(fs/2),'bandpass');

figure('Name','precision spot-check')
for j = 1:nShow
    k = showIdx(j); idx = ripplesClean(k,2)+sRelFull;
    if idx(1)<1 || idx(end)>nSamplesTotal, continue; end
    xx=double(m.Data.x(chOrder,idx))*lsb; car=median(xx,1);
    pyr=xx(find(chOrder==pyrCh),:)-car; rad=xx(find(chOrder==radCh),:)-car;
    fp = filtfilt(bR,aR,pyr); fq = filtfilt(bS,aS,rad);   % filter full, then crop
    subplot(ceil(nShow/2),2,j)
    plot(t, fp(crop),'k'); hold on; plot(t, fq(crop),'b');
    xlim(plotWin); ylabel(['#' num2str(k)])
    if j==1, title('ripple(k) + SW(b)'); end
end


%% ===================================================================
%  G. SWR-confirming figures + feature characterization
%  Run after the pipeline (needs ripplesClean, m, fs, chOrder, pyrCh,
%  radCh, rippleBand, swBand, nSamplesTotal). Each panel makes a specific
%  argument that the curated events are genuine CA1 sharp-wave ripples.
%% ===================================================================


%% G1. Mean ripple waveform (trough-aligned) + mean sharp wave ---------
%  Aligning each event to its central ripple trough lets the cycles add
%  coherently (envelope-peak alignment averages the oscillation to ~0).
%  A clean, symmetric, damped oscillation = a real rhythmic ripple.
wfWin   = [-0.05 0.05];
outRel  = round(wfWin(1)*fs):round(wfWin(2)*fs);   % output window, rel. to trough
tWf     = outRel/fs;
alignWin= round(0.008*fs);                         % +/-8 ms search for central trough
padRel  = round((wfWin(1)-0.02)*fs):round((wfWin(2)+0.02)*fs);  % padded read window
c0      = find(padRel==0);                         % index of t=0 in padded window

[bR,aR] = butter(3, rippleBand/(fs/2),'bandpass');
[bS,aS] = butter(3, swBand/(fs/2),'bandpass');
chR = unique([pyrCh radCh hpcCh]);
iP = find(chR==pyrCh); iS = find(chR==radCh); iRef = find(ismember(chR,hpcCh));

ev = ripplesClean(:,2);
Wr = []; Wraw = []; Wsw = [];
for k = 1:numel(ev)
    idx = ev(k)+padRel;
    if idx(1)<1 || idx(end)>nSamplesTotal, continue; end
    xx  = double(m.Data.x(chR,idx))*lsb;
    car = median(xx(iRef,:),1);
    pyr = xx(iP,:)-car; rad = xx(iS,:)-car;
    fr  = filtfilt(bR,aR,pyr);
    frsw= filtfilt(bS,aS,rad);

    seg = (c0-alignWin):(c0+alignWin);
    [~,mi] = min(fr(seg)); tr = seg(mi);           % trough index in padded window
    oi = tr + outRel;                              % extract output window around trough
    if oi(1)<1 || oi(end)>numel(padRel), continue; end

    Wr(end+1,:)   = fr(oi);     %#ok<SAGROW>
    Wraw(end+1,:) = pyr(oi);    %#ok<SAGROW>
    Wsw(end+1,:)  = frsw(oi);   %#ok<SAGROW>
end

mR = mean(Wr,1); seR = std(Wr,[],1)/sqrt(size(Wr,1));
figure('Name','G1 mean ripple waveform')
subplot(1,2,1)
fill([tWf fliplr(tWf)]*1e3, [mR+seR fliplr(mR-seR)], [.8 .8 1],'EdgeColor','none'); hold on
plot(tWf*1e3, mR, 'b')
xlabel('time from ripple trough (ms)'); ylabel('\muV')
title(sprintf('mean ripple (trough-aligned), n=%d', size(Wr,1)))
subplot(1,2,2)
plot(tWf*1e3, mean(Wraw,1), 'k'); hold on; plot(tWf*1e3, mean(Wsw,1), 'b')
legend('pyr broadband','rad sharp wave'); xlabel('time from ripple trough (ms)')
ylabel('\muV'); title('broadband + sharp wave')


%% G2. Event-averaged spectrogram + population PSD --------------------
%  THE key argument vs gamma/MUA: a transient, NARROW-band peak at the
%  ripple frequency riding on the 1/f background, averaged over all events.
tfWin = [-0.15 0.15];
sRel  = round(tfWin(1)*fs):round(tfWin(2)*fs);
chR = unique([pyrCh hpcCh]); iP = find(chR==pyrCh); iRef = find(ismember(chR,hpcCh));
ev  = ripplesClean(:,2);

nwin = round(0.025*fs); nov = nwin - round(0.002*fs); nf = 0:4:300;
S = 0; cnt = 0;
for k = 1:numel(ev)
    idx = ev(k)+sRel;
    if idx(1)<1 || idx(end)>nSamplesTotal, continue; end
    xx = double(m.Data.x(chR,idx))*lsb; pyr = xx(iP,:)-median(xx(iRef,:),1);
    [s,fS,tS] = spectrogram(pyr, nwin, nov, nf, fs);
    P = abs(s).^2;
    if cnt==0, S = P; else, S = S + P; end
    cnt = cnt + 1;
end
S = S/cnt;

% population PSD: event window vs flank baseline, averaged over events
pWin = [-0.3 0.3]; sRelP = round(pWin(1)*fs):round(pWin(2)*fs); ttp = sRelP/fs;
evM = ttp>=-0.04 & ttp<=0.04; baM = ttp>=-0.30 & ttp<=-0.15;
nfft = 2^nextpow2(round(0.08*fs)); Pe = 0; Pb = 0; cnt2 = 0;
for k = 1:numel(ev)
    idx = ev(k)+sRelP;
    if idx(1)<1 || idx(end)>nSamplesTotal, continue; end
    xx = double(m.Data.x(chR,idx))*lsb; pyr = xx(iP,:)-median(xx(iRef,:),1);
    [pe,fp] = pwelch(pyr(evM),[],[],nfft,fs);
    [pb,~]  = pwelch(pyr(baM),[],[],nfft,fs);
    Pe = Pe + pe; Pb = Pb + pb; cnt2 = cnt2 + 1;
end
Pe = Pe/cnt2; Pb = Pb/cnt2;

tAx   = tS + tfWin(1);
Sbase = mean(S(:, tAx < -0.05), 2);            % pre-event baseline per frequency
Snorm = 10*log10(S ./ Sbase);                  % dB relative to baseline -> transient pops
figure('Name','G2 event TF + population PSD')
subplot(1,2,1)
imagesc(tAx, fS, Snorm); axis xy; ylim([40 280]); caxis([-2 8])
xlabel('time from peak (s)'); ylabel('Hz'); colorbar; colormap('viridis')
title(sprintf('mean spectrogram (dB vs baseline), n=%d', cnt))
subplot(1,2,2)
plot(fp, 10*log10(Pe), 'k'); hold on; plot(fp, 10*log10(Pb), 'Color',[.6 .6 .6])
xlim([0 300]); xlabel('Hz'); ylabel('dB'); legend('event','baseline')
title('population PSD (narrow peak on 1/f)')


%% G3. Laminar double dissociation -----------------------------------
%  Ripple power should peak in the PYRAMIDAL layer; the sharp-wave sink in
%  RADIATUM. Two distinct generators in two layers = CA1 SWR signature.
nC = numel(chOrder);
iPyr = find(chOrder==pyrCh); iRad = find(chOrder==radCh);
lWin = [-0.05 0.05]; sRel = round(lWin(1)*fs):round(lWin(2)*fs); tL = sRel/fs;
[bR,aR] = butter(3, rippleBand/(fs/2),'bandpass');
[bS,aS] = butter(3, swBand/(fs/2),'bandpass');
ev = ripplesClean(:,2);

ripEnvDepth = zeros(nC,1); swDepth = zeros(nC,numel(sRel)); cnt = 0;
for k = 1:numel(ev)
    idx = ev(k)+sRel;
    if idx(1)<1 || idx(end)>nSamplesTotal, continue; end
    xx = double(m.Data.x(chOrder,idx))*lsb; xx = xx - median(xx,1);
    fr = filtfilt(bR,aR,xx')'; env = abs(hilbert(fr')');   % env: nC x time
    ripEnvDepth = ripEnvDepth + mean(env,2);
    swDepth = swDepth + filtfilt(bS,aS,xx')';
    cnt = cnt + 1;
end
ripEnvDepth = ripEnvDepth/cnt; swDepth = swDepth/cnt;

rpN = ripEnvDepth / max(ripEnvDepth);
swN = -min(swDepth(:,tL>-0.02 & tL<0.02),[],2); swN = swN/max(swN);

figure('Name','G3 laminar dissociation')
subplot(1,3,1)
plot(ripEnvDepth, 1:nC, 'k.-'); axis ij; hold on; plot(ripEnvDepth(iPyr), iPyr, 'ro')
xlabel('mean ripple envelope'); ylabel('depth rank'); title('ripple power -> pyr')
subplot(1,3,2)
imagesc(tL, 1:nC, swDepth); axis xy; hold on; plot(0, iRad, 'rx')
xlabel('time (s)'); title('sharp wave (depth)'); colormap('viridis'); colorbar
subplot(1,3,3)
plot(rpN, 1:nC, 'k'); axis ij; hold on; plot(swN, 1:nC, 'b')
plot(xlim, [iPyr iPyr], 'r:')
legend('ripple env','SW sink','pyr'); xlabel('normalized'); ylabel('depth rank')
title('lamination')


%% G4. Temporal structure: raster, rate, IEI, autocorrelogram ---------
[rt,ord] = sort(rippleTimesClean(:));
evS = ripplesClean(ord,2);
durTot = nSamplesTotal/fs;
iPyrG = find(chOrder==pyrCh);

% raster ACROSS DEPTH: per-event ripple-power profile binned in real time.
% Events appear as vertical streaks confined to the pyramidal layer.
w = round(0.02*fs); sRelR = -w:w;
[bRf,aRf] = butter(3, rippleBand/(fs/2),'bandpass');
binS = 2; nCol = ceil(durTot/binS);
Mtime = zeros(nC,nCol);
for k = 1:numel(evS)
    idx = evS(k)+sRelR;
    if idx(1)<1 || idx(end)>nSamplesTotal, continue; end
    xx = double(m.Data.x(chOrder,idx))*lsb; xx = xx - median(xx,1);
    fr = filtfilt(bRf,aRf,xx')'; env = abs(hilbert(fr')');
    col = max(1, ceil(rt(k)/binS));
    Mtime(:,col) = Mtime(:,col) + max(env,[],2);
end

figure('Name','G4 temporal structure')
subplot(4,1,1)
imagesc((1:nCol)*binS/60, 1:nC, Mtime); axis xy; colormap('viridis')
hold on; plot(xlim, [iPyrG iPyrG], 'r:')
xlabel('time (min)'); ylabel('channel (depth)')
title(sprintf('event raster across depth (%d s bins), n=%d', binS, numel(evS)))

subplot(4,1,2)
bw = 30; edges = 0:bw:durTot; rate = histcounts(rt, edges)/bw;
stairs(edges(1:end-1)/60, rate, 'k'); xlim([0 durTot/60])
xlabel('time (min)'); ylabel('events/s'); title(sprintf('rate (%d s bins)', bw))

subplot(4,1,3)
iei = diff(rt);
histogram(iei(iei<3), 0:0.025:3); xlabel('inter-event interval (s)'); ylabel('count')
title('IEI (refractoriness + typical spacing)')

subplot(4,1,4)
maxLag = 2; binAC = 0.02; d = [];
for i = 1:numel(rt)
    dd = rt - rt(i); dd = dd(abs(dd)<=maxLag & dd~=0); d = [d; dd]; %#ok<AGROW>
end
edgesAC = -maxLag:binAC:maxLag; ac = histcounts(d, edgesAC); ctr = edgesAC(1:end-1)+binAC/2;
bar(ctr, ac, 'k'); xlabel('lag (s)'); ylabel('count'); title('event autocorrelogram')


%% G5. Sharp-wave / ripple coupling + intra-ripple frequency (chirp) ---
%  (1) ripple amplitude vs sharp-wave amplitude  -> coupling strength
%  (2) ripple-peak minus SW-trough lag           -> timing
%  (3) instantaneous ripple frequency vs time     -> chirp (deceleration)
cWin = [-0.1 0.1]; sRel = round(cWin(1)*fs):round(cWin(2)*fs);
[bR,aR] = butter(3, rippleBand/(fs/2),'bandpass');
[bS,aS] = butter(3, swBand/(fs/2),'bandpass');
[bW,aW] = butter(3, [120 280]/(fs/2),'bandpass');   % wider band so the chirp is not squashed
chR = unique([pyrCh radCh hpcCh]); iP=find(chR==pyrCh); iS=find(chR==radCh); iRef=find(ismember(chR,hpcCh));
ev = ripplesClean(:,2); N = numel(ev);
ripAmp=nan(N,1); swAmp=nan(N,1); lag=nan(N,1);
ifWin=[-0.025 0.025]; sIf=round(ifWin(1)*fs):round(ifWin(2)*fs); tIf=sIf/fs;
IF = nan(N,numel(sIf));
for k = 1:N
    idx = ev(k)+sRel; if idx(1)<1 || idx(end)>nSamplesTotal, continue; end
    xx = double(m.Data.x(chR,idx))*lsb; car = median(xx(iRef,:),1);
    pyr = xx(iP,:)-car; rad = xx(iS,:)-car;
    fr = filtfilt(bR,aR,pyr); env = abs(hilbert(fr));
    sw = filtfilt(bS,aS,rad);
    [ripAmp(k),ip] = max(env);               % ripple peak envelope (uV) + its index
    [swMin,is]     = min(sw);                 % SW trough (uV) + its index
    swAmp(k) = -swMin;                        % sink depth as a positive amplitude
    lag(k)   = (ip-is)/fs;                    % ripple peak - SW trough (s)
    frW = filtfilt(bW,aW,pyr); phW = unwrap(angle(hilbert(frW)));
    instf = diff(phW)/(2*pi)*fs;              % instantaneous frequency (length N-1)
    oi = ip + sIf;
    if oi(1)>=1 && oi(end)<=numel(instf), IF(k,:) = instf(oi); end
end
v = ~isnan(ripAmp) & ~isnan(swAmp); R = corrcoef(swAmp(v),ripAmp(v)); rSW = R(1,2);
mIF = mean(IF,1,'omitnan'); sdIF = std(IF,[],1,'omitnan'); nIF = sum(~isnan(IF),1); seIF = sdIF./sqrt(nIF);

figure('Name','G5 SW-ripple coupling + chirp')
subplot(1,3,1)
plot(swAmp, ripAmp, '.', 'Color',[.3 .3 .8])
xlabel('sharp-wave amplitude (\muV)'); ylabel('ripple amplitude (\muV)')
title(sprintf('coupling: r = %.2f (n=%d)', rSW, sum(v)))
subplot(1,3,2)
histogram(lag*1e3, -50:5:50); xlabel('ripple peak - SW trough (ms)'); ylabel('count')
title(sprintf('SW-ripple lag (median %.0f ms)', median(lag*1e3,'omitnan')))
subplot(1,3,3)
fill([tIf fliplr(tIf)]*1e3, [mIF+seIF fliplr(mIF-seIF)], [.8 .8 1], 'EdgeColor','none','HandleVisibility','off'); hold on
plot(tIf*1e3, mIF, 'b'); xlabel('time from ripple peak (ms)'); ylabel('inst. frequency (Hz)')
title('intra-ripple frequency (chirp)')


%% G6. Singlets / duplets / triplets (ripple bursts) ------------------
%  Group ripples whose peaks fall within maxGap into a burst, then count
%  burst sizes. maxGap is the key parameter - tune it to your IEI histogram.
maxGap = 0.12;   % s
[rt,ord] = sort(rippleTimesClean(:));
g = [Inf; diff(rt)];
clustID   = cumsum(g > maxGap);            % burst index per event (time order)
clustSize = accumarray(clustID, 1);        % # ripples in each burst
bsSorted  = clustSize(clustID);
burstSize_clean = zeros(size(rt));
burstSize_clean(ord) = bsSorted;           % burst size aligned to ripplesClean rows (for downstream)

fprintf('Bursts (gap<%.0f ms): %d singlets, %d duplets, %d triplets, %d quads+\n', ...
    1e3*maxGap, sum(clustSize==1), sum(clustSize==2), sum(clustSize==3), sum(clustSize>=4));

[bR,aR] = butter(3, rippleBand/(fs/2),'bandpass');
[bS,aS] = butter(3, swBand/(fs/2),'bandpass');
chR = unique([pyrCh radCh hpcCh]); iP=find(chR==pyrCh); iS=find(chR==radCh); iRef=find(ismember(chR,hpcCh));

figure('Name','G6 duplets/triplets')
subplot(1,3,1)
histogram(clustSize, 0.5:1:max(6,max(clustSize))+0.5)
xlabel('ripples per burst'); ylabel('# bursts'); title(sprintf('burst size (gap<%d ms)', 1e3*maxGap))

exTargets = [2 3];
for ti = 1:numel(exTargets)
    subplot(1,3,1+ti)
    ci = find(clustSize==exTargets(ti), 1, 'first');
    if isempty(ci), title(sprintf('no %d-ripple burst', exTargets(ti))); continue; end
    ts = rt(clustID==ci); tc = mean(ts);
    sR = round((tc-0.15)*fs):round((tc+0.15)*fs);
    if sR(1)<1 || sR(end)>nSamplesTotal, continue; end
    xx = double(m.Data.x(chR,sR))*lsb; car = median(xx(iRef,:),1);
    pyr = xx(iP,:)-car; rad = xx(iS,:)-car; tb = sR/fs - tc;
    plot(tb, filtfilt(bR,aR,pyr), 'k'); hold on; plot(tb, filtfilt(bS,aS,rad), 'b')
    yl = ylim; plot(ts-tc, yl(2)*ones(size(ts)), 'rv', 'MarkerFaceColor','r')
    xlabel('time (s)'); title(sprintf('example %d-ripple burst', exTargets(ti)))
end


%% G5b. SW-ripple amplitude coupling on UNCURATED detections ----------
%  The curated coupling (G5) is range-restricted by the swZ>2 gate. Recompute
%  on ALL detected events for an unbiased estimate (curated highlighted).
cWin = [-0.1 0.1]; sRel = round(cWin(1)*fs):round(cWin(2)*fs);
[bR,aR] = butter(3, rippleBand/(fs/2),'bandpass');
[bS,aS] = butter(3, swBand/(fs/2),'bandpass');
chR = unique([pyrCh radCh hpcCh]); iP=find(chR==pyrCh); iS=find(chR==radCh); iRef=find(ismember(chR,hpcCh));
evA = ripples(:,2); NA = numel(evA);
ripAmpA = nan(NA,1); swAmpA = nan(NA,1);
for k = 1:NA
    idx = evA(k)+sRel; if idx(1)<1 || idx(end)>nSamplesTotal, continue; end
    xx = double(m.Data.x(chR,idx))*lsb; car = median(xx(iRef,:),1);
    pyr = xx(iP,:)-car; rad = xx(iS,:)-car;
    ripAmpA(k) = max(abs(hilbert(filtfilt(bR,aR,pyr))));
    swAmpA(k)  = -min(filtfilt(bS,aS,rad));
end
isC = ismember(evA, ripplesClean(:,2));
vA = ~isnan(ripAmpA) & ~isnan(swAmpA);
RA = corrcoef(swAmpA(vA), ripAmpA(vA)); rA = RA(1,2);
vC = vA & isC; RC = corrcoef(swAmpA(vC), ripAmpA(vC)); rC = RC(1,2);

figure('Name','G5b coupling (uncurated)')
plot(swAmpA(vA & ~isC), ripAmpA(vA & ~isC), '.', 'Color',[.65 .65 .65]); hold on
plot(swAmpA(vC), ripAmpA(vC), '.', 'Color',[.85 .33 .1])
xlabel('sharp-wave amplitude (\muV)'); ylabel('ripple amplitude (\muV)')
legend('all detected','curated')
title(sprintf('coupling: all r=%.2f (n=%d) | curated r=%.2f (n=%d)', rA, sum(vA), rC, sum(vC)))


%% G6b. Burst maxGap sweep + sustained-vs-separate sharp wave ----------
[rt,ord] = sort(rippleTimesClean(:));
gaps = diff(rt);

% (1) sweep: how burst composition depends on the (arbitrary) gap threshold
gapVec = 0.05:0.025:0.40;
fracSing=zeros(size(gapVec)); fracDup=fracSing; fracTrip=fracSing;
for gi = 1:numel(gapVec)
    cID = cumsum([true; gaps>gapVec(gi)]);
    cS  = accumarray(cID,1);
    fracSing(gi)=mean(cS==1); fracDup(gi)=mean(cS==2); fracTrip(gi)=mean(cS>=3);
end

% (2) sub-classify multi-ripple bursts: does the radiatum sink stay elevated
%     BETWEEN consecutive ripple peaks (one prolonged SW) or return toward
%     baseline (separate SW-R complexes)?
maxGap = 0.12;
cID = cumsum([true; gaps>maxGap]); cS = accumarray(cID,1);
[bSf,aSf] = butter(3, swBand/(fs/2),'bandpass');
chR = unique([radCh hpcCh]); iS = find(chR==radCh); iRef = find(ismember(chR,hpcCh));
multi = find(cS>=2); nSust = 0; nSep = 0;
for b = 1:numel(multi)
    pk = rt(cID==multi(b)); sepPairs = 0;
    for j = 1:numel(pk)-1
        t1=pk(j); t2=pk(j+1); tm=(t1+t2)/2;
        sR = round((t1-0.03)*fs):round((t2+0.03)*fs);
        if sR(1)<1 || sR(end)>nSamplesTotal, continue; end
        xx = double(m.Data.x(chR,sR))*lsb; car = median(xx(iRef,:),1);
        sw = filtfilt(bSf,aSf, xx(iS,:)-car); tax = sR/fs;
        sw1=sw(find(tax>=t1,1)); sw2=sw(find(tax>=t2,1)); swm=sw(find(tax>=tm,1));
        depth = min(sw1,sw2);                       % deeper of the two troughs
        if ~(depth<0 && swm <= 0.5*depth), sepPairs = sepPairs+1; end   % mid returned to baseline
    end
    if sepPairs==0, nSust = nSust+1; else nSep = nSep+1; end
end

figure('Name','G6b burst sweep + sub-type')
subplot(1,2,1)
plot(gapVec*1e3,fracSing,'k', gapVec*1e3,fracDup,'b', gapVec*1e3,fracTrip,'r')
xlabel('maxGap (ms)'); ylabel('fraction of bursts'); legend('singlet','duplet','triplet+')
title('burst composition vs maxGap')
subplot(1,2,2)
bar([nSust nSep]); set(gca,'xticklabel',{'sustained (1 SW)','separate SWs'})
ylabel('# multi-ripple bursts'); title(sprintf('sub-type (gap<%d ms, n=%d)', 1e3*maxGap, numel(multi)))


%% G7. Example single SWRs: ripple (pyr) riding on the sharp wave (rad)
%  Per event, three large traces: the slow sharp wave from str. radiatum
%  (blue), the band-passed ripple from the pyramidal layer (red), and their
%  SUM (black) = the ripple visibly sitting on the sharp-wave deflection.
%  (Single events: the mean cancels the ripple. cf. Buzsaki 1992; Liu 2022.)
nEx = 3;
exWin = [-0.06 0.06]; sRel = round(exWin(1)*fs):round(exWin(2)*fs); tM = sRel/fs;
[bR,aR] = butter(3, rippleBand/(fs/2),'bandpass');
[bS,aS] = butter(3, swBand/(fs/2),'bandpass');
chR = unique([pyrCh radCh hpcCh]); iP=find(chR==pyrCh); iS=find(chR==radCh); iRef=find(ismember(chR,hpcCh));
[~,ordZ] = sort(ripplesClean(:,4),'descend');
exIdx = ordZ(round(linspace(1, min(40,size(ripplesClean,1)), nEx)));
exEv = ripplesClean(exIdx,2);

figure('Name','G7 example SWRs: ripple on sharp wave')
for e = 1:numel(exEv)
    idx = exEv(e)+sRel; if idx(1)<1 || idx(end)>nSamplesTotal, continue; end
    xx = double(m.Data.x(chR,idx))*lsb; car = median(xx(iRef,:),1);
    pyrRaw = xx(iP,:)-car; radRaw = xx(iS,:)-car;
    pyrRip = filtfilt(bR,aR,pyrRaw);
    radSW  = filtfilt(bS,aS,radRaw);
    comp   = radSW + pyrRip;                       % ripple riding on the sharp wave
    gap = 1.5*max(abs(radSW));
    subplot(1,numel(exEv),e); hold on
    plot(tM*1e3, radSW,          'b')
    plot(tM*1e3, pyrRip + gap,   'r')
    plot(tM*1e3, comp   + 2*gap, 'k')
    plot([0 0], [min(radSW)-0.5*gap, 2.5*gap+max(abs(comp))], 'k:')
    set(gca,'ytick',[0 gap 2*gap],'yticklabel',{'sharp wave (rad)','ripple (pyr)','sum'})
    xlim(exWin*1e3); xlabel('time from ripple peak (ms)')
    title(sprintf('event z = %.1f', ripplesClean(exIdx(e),4)))
end
