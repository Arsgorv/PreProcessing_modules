function R = detect_ripples(datapath, varargin)
%DETECT_RIPPLES  Detect hippocampal sharp-wave ripples from one session.
%
%   R = detect_ripples(datapath)
%   R = detect_ripples(datapath, 'Name',Value, ...)
%
% Expects datapath to contain 'continuous.dat' and 'continuous.session.mat'
% (CellExplorer-style; fs/nChannels/lsb read from session.extracellular).
%
% Pipeline (matches the validated exploration script): depth-order channels,
% localize CA1 pyramidal (max ripple power) and radiatum (sharp-wave sink),
% pick the ripple band from the data, detect (chunked, MAD-z dual threshold,
% merged), require a co-occurring sharp wave, tag each event by epoch, and
% save <outdir>/<session>_ripples.mat (+ a QC figure). Any auto step can be
% overridden. ALWAYS check the QC figure per session.
%
% Key options (defaults):
%   outdir (=datapath), hpcCh (=all), badCh ([]), probeType ('NP1.0'),
%   pyrCh ([]=auto), radCh ([]=auto), rippleBand ([]=auto),
%   thrHigh (5), thrLow (2.5), minDur (.02), maxDur (.2), minInterval (.03),
%   swBand ([2 40]), swZthr (2), saveFig (true)
%   epochs ([] or [nEpoch x 2] start/stop in SECONDS) -- e.g. NREM intervals
%   epochNames ({} cellstr, one per epoch row)
%   restrict (false) -- if true, keep only ripples falling inside 'epochs'
%
% Output struct R has peakTime, startStop, duration, events, per-event
% features, epoch label per ripple, per-epoch rate, and the params used.

p = inputParser;
p.addParameter('outdir', datapath);
p.addParameter('hpcCh', []);     p.addParameter('badCh', []);
p.addParameter('probeType', 'NP1.0');
p.addParameter('pyrCh', []);     p.addParameter('radCh', []);
p.addParameter('rippleBand', []);
p.addParameter('thrHigh', 5);    p.addParameter('thrLow', 2.5);
p.addParameter('minDur', 0.020); p.addParameter('maxDur', 0.200);
p.addParameter('minInterval', 0.030);
p.addParameter('swBand', [2 40]); p.addParameter('swZthr', 2);
p.addParameter('epochs', []);    p.addParameter('epochNames', {});
p.addParameter('restrict', false);
p.addParameter('tsdUnits', 1e4); % intervalSet/ts time units: time_sec * tsdUnits
p.addParameter('saveFig', true);
p.parse(varargin{:}); o = p.Results;

% ---- session ----
datFile = fullfile(datapath,'continuous.dat');
[fs, nCh, lsb, yc0] = load_session_info(datapath);   % session.mat or structure.oebin
fi  = dir(datFile); nST = fi.bytes/2/nCh;
m   = memmapfile(datFile,'Format',{'int16',[nCh nST],'x'});
parts = regexp(datapath,'[\\/]','split'); parts = parts(~cellfun(@isempty,parts));
sessionName = parts{end};

if isempty(o.hpcCh), o.hpcCh = 1:nCh; end
hpcCh = setdiff(o.hpcCh, o.badCh);

% ---- depth order ----
if ~isempty(yc0)
    yc = yc0;
else
    if strcmpi(o.probeType,'NP2.0'), vp = 15; else, vp = 20; end
    yc = vp*floor((0:nCh-1)'/2);
end
[~,dsort] = sort(yc(hpcCh)); chOrder = hpcCh(dsort); nC = numel(chOrder);

% ---- pyramidal channel (max ripple-band power over sample windows) ----
pyrCh = o.pyrCh;
if isempty(pyrCh)
    segLen = round(20*fs); nSamp = 8;
    starts = round(linspace(segLen, max(segLen+1, nST-2*segLen), nSamp));
    [bL,aL] = butter(3,[90 200]/(fs/2),'bandpass'); pw = zeros(nC,1);
    for s = 1:nSamp
        idx = starts(s):starts(s)+segLen-1; idx = idx(idx>=1 & idx<=nST);
        xx = double(m.Data.x(chOrder,idx))*lsb; xx = xx - median(xx,1);
        fr = filtfilt(bL,aL,xx')'; env = abs(hilbert(fr')');
        pw = pw + mean(env.^2,2);
    end
    [~,iP] = max(pw); pyrCh = chOrder(iP);
end
iPyr = find(chOrder==pyrCh);

% ---- ripple band (PSD peak from a provisional detection) ----
rippleBand = o.rippleBand;
if isempty(rippleBand)
    prov = local_detect(m,nST,fs,lsb,chOrder,pyrCh,hpcCh,[90 200], ...
        o.thrHigh,o.thrLow,o.minDur,o.maxDur,o.minInterval);
    rippleBand = local_band(m,nST,fs,lsb,chOrder,pyrCh,hpcCh,prov(:,2));
end

% ---- radiatum channel (sharp-wave sink from provisional events) ----
radCh = o.radCh;
if isempty(radCh)
    if ~exist('prov','var')
        prov = local_detect(m,nST,fs,lsb,chOrder,pyrCh,hpcCh,rippleBand, ...
            o.thrHigh,o.thrLow,o.minDur,o.maxDur,o.minInterval);
    end
    radCh = local_radch(m,nST,fs,lsb,chOrder,hpcCh,o.swBand,prov(:,2));
end
iRad = find(chOrder==radCh);

% ---- final detection ----
ripples = local_detect(m,nST,fs,lsb,chOrder,pyrCh,hpcCh,rippleBand, ...
    o.thrHigh,o.thrLow,o.minDur,o.maxDur,o.minInterval);

% ---- curation: features + sharp-wave gate ----
[peakZ,swZ,fPeak,nCyc,ratioHi] = local_features(m,nST,fs,lsb,chOrder, ...
    pyrCh,radCh,hpcCh,rippleBand,o.swBand,ripples);
good = peakZ>o.thrHigh & swZ>o.swZthr & nCyc>=4 & ratioHi>2;
ripplesClean = ripples(good,:);
pkT = ripplesClean(:,2)/fs;

% ---- epochs: tag each ripple, per-epoch rate, optional restriction ----
ep = o.epochs;
if ~isempty(ep) && exist('intervalSet')~=0 && isa(ep,'intervalSet')
    ep = [Start(ep) End(ep)] / o.tsdUnits;     % intervalSet -> [start stop] seconds
end
epIdx = zeros(numel(pkT),1);
epRate = []; epDur = [];
if ~isempty(ep)
    epIdx = local_epochidx(pkT, ep);
    epDur = ep(:,2)-ep(:,1);
    epRate = accumarray([epIdx(epIdx>0); size(ep,1)], [ones(sum(epIdx>0),1); 0], [size(ep,1) 1]) ./ max(epDur,eps);
    if o.restrict
        keep = epIdx>0;
        ripplesClean = ripplesClean(keep,:); pkT = pkT(keep);
        peakZ_g = peakZ(good); swZ_g = swZ(good); fPeak_g = fPeak(good); nCyc_g = nCyc(good);
        peakZ_g = peakZ_g(keep); swZ_g = swZ_g(keep); fPeak_g = fPeak_g(keep); nCyc_g = nCyc_g(keep);
        epIdx = epIdx(keep);
    end
end
if ~exist('peakZ_g','var')
    peakZ_g = peakZ(good); swZ_g = swZ(good); fPeak_g = fPeak(good); nCyc_g = nCyc(good);
end

% ---- output ----
R = struct();
R.sessionName = sessionName; R.datapath = datapath; R.fs = fs;
R.pyrCh = pyrCh; R.radCh = radCh; R.rippleBand = rippleBand;
R.peakTime  = pkT;
R.startStop = [ripplesClean(:,1) ripplesClean(:,3)]/fs;
R.duration  = ripplesClean(:,5);
R.events    = ripplesClean;          % [startSample peakSample stopSample peakZ dur]
R.features  = struct('peakZ',peakZ_g,'swZ',swZ_g,'peakFreq',fPeak_g,'nCycles',nCyc_g);
R.epoch     = epIdx;                  % per-ripple epoch index (0 = outside all epochs)
R.epochs    = ep; R.epochNames = o.epochNames;
R.epochRate = epRate; R.epochDur = epDur;
R.nDetected = size(ripples,1); R.nCurated = numel(pkT);
R.recDurSec = nST/fs;
R.params    = o;

% intervalSet / ts for TStoolbox restriction workflows, if available:
%   Restrict(anyTSD, R.rippleIS)  -> that signal during ripples
%   Restrict(R.rippleTS, stateIS) -> ripple peaks in a given state
if exist('intervalSet') ~= 0
    try
        uu = o.tsdUnits;
        R.rippleIS = intervalSet(R.startStop(:,1)*uu, R.startStop(:,2)*uu);
        R.rippleTS = ts(R.peakTime*uu);
        R.tsdUnits = uu;
    catch ME
        warning('detect_ripples:tsd','Could not build intervalSet/ts (%s).', ME.message);
    end
end

if ~exist(o.outdir,'dir'), mkdir(o.outdir); end
save(fullfile(o.outdir,[sessionName '_ripples.mat']),'R','-v7.3');
fprintf('%s: pyr %d, rad %d, band [%d %d] Hz | %d/%d curated | %.3f/s (whole rec)\n', ...
    sessionName, pyrCh, radCh, rippleBand(1), rippleBand(2), R.nCurated, R.nDetected, R.nCurated/(nST/fs));
if ~isempty(epRate)
    for e = 1:numel(epRate)
        nm = sprintf('epoch %d',e); if numel(o.epochNames)>=e, nm = o.epochNames{e}; end
        fprintf('   %-14s rate %.3f /s  (%.1f s)\n', nm, epRate(e), epDur(e));
    end
end

if o.saveFig
    local_qcfig(m,nST,fs,lsb,chOrder,iPyr,iRad,pyrCh,radCh,rippleBand,o.swBand, ...
        ripplesClean(:,2), R, fullfile(o.outdir,[sessionName '_ripples_QC.png']));
end
end


% ===================================================================
function ripples = local_detect(m,nST,fs,lsb,chOrder,detectCh,refList,band, ...
        thrHigh,thrLow,minDur,maxDur,minInterval)
[b,a] = butter(3, band/(fs/2),'bandpass');
chans = unique([detectCh refList]); di = find(chans==detectCh); ri = ismember(chans,refList);
chunk = round(60*fs); ov = round(1*fs); ripples = []; cs = 1;
while cs <= nST
    ce = min(nST,cs+chunk-1); rs = max(1,cs-ov); re = min(nST,ce+ov);
    xc = double(m.Data.x(chans,rs:re))*lsb;
    x  = xc(di,:) - median(xc(ri,:),1);
    xf = filtfilt(b,a,x); env = abs(hilbert(xf));
    mu = median(env); sd = median(abs(env-mu))*1.4826;
    if sd==0 || isnan(sd), cs = ce+1; continue; end
    z = (env-mu)/sd; aH = z>thrHigh; aL = z>thrLow;
    d = diff([0 aH 0]); ps = find(d==1); pe = find(d==-1)-1;
    ks = cs-rs+1; ke = ce-rs+1;
    for i = 1:numel(ps)
        s1 = ps(i); while s1>1 && aL(s1), s1 = s1-1; end
        s2 = pe(i); while s2<numel(aL) && aL(s2), s2 = s2+1; end
        [pv,im] = max(z(s1:s2)); pk = s1+im-1;
        if pk<ks || pk>ke, continue; end
        dur = (s2-s1+1)/fs; if dur<minDur || dur>maxDur, continue; end
        ripples = [ripples; rs+s1-1, rs+pk-1, rs+s2-1, pv, dur]; %#ok<AGROW>
    end
    cs = ce+1;
end
if ~isempty(ripples)
    ripples = sortrows(ripples,1); mg = ripples(1,:);
    for i = 2:size(ripples,1)
        if (ripples(i,1)-mg(end,3))/fs < minInterval
            mg(end,3) = ripples(i,3);
            if ripples(i,4) > mg(end,4), mg(end,2) = ripples(i,2); mg(end,4) = ripples(i,4); end
            mg(end,5) = (mg(end,3)-mg(end,1))/fs;
        else
            mg = [mg; ripples(i,:)]; %#ok<AGROW>
        end
    end
    ripples = mg;
end
end


% ===================================================================
function band = local_band(m,nST,fs,lsb,chOrder,pyrCh,refList,ev)
chans = unique([pyrCh refList]); pidx = find(chans==pyrCh); ri = ismember(chans,refList);
w = [-0.3 0.3]; sRel = round(w(1)*fs):round(w(2)*fs); tt = sRel/fs;
em = tt>=-0.04 & tt<=0.04; bm = tt>=-0.30 & tt<=-0.15;
nfft = 2^nextpow2(round(0.08*fs)); Pe = 0; Pb = 0; n = 0;
if numel(ev)>500, ev = ev(randperm(numel(ev),500)); end
for k = 1:numel(ev)
    idx = ev(k)+sRel; if idx(1)<1 || idx(end)>nST, continue; end
    xx = double(m.Data.x(chans,idx))*lsb; py = xx(pidx,:) - median(xx(ri,:),1);
    [pe,f] = pwelch(py(em),[],[],nfft,fs); [pb,~] = pwelch(py(bm),[],[],nfft,fs);
    Pe = Pe+pe; Pb = Pb+pb; n = n+1;
end
ratio = (Pe/n)./((Pb/n)+eps);
sel = f>=100 & f<=250; fa = f(sel); rr = ratio(sel);
isp = false(numel(rr),1); isp(2:end-1) = rr(2:end-1)>rr(1:end-2) & rr(2:end-1)>rr(3:end);
loc = find(isp); if isempty(loc), [~,ip] = max(rr); else, [~,j] = max(rr(loc)); ip = loc(j); end
fp = fa(ip); nyq = fs/2;
band = round([0.80*fp 1.15*fp]/5)*5; band(2) = min(band(2), floor(0.9*nyq));
end


% ===================================================================
function radCh = local_radch(m,nST,fs,lsb,chOrder,refList,swBand,ev)
nC = numel(chOrder); [bS,aS] = butter(3, swBand/(fs/2),'bandpass');
chans = unique([chOrder(:)' refList]); [~,ci] = ismember(chOrder,chans); ri = ismember(chans,refList);
w = [-0.05 0.05]; sRel = round(w(1)*fs):round(w(2)*fs); tt = sRel/fs; cm = tt>-0.015 & tt<0.015;
prof = zeros(nC,numel(sRel)); n = 0;
if numel(ev)>500, ev = ev(randperm(numel(ev),500)); end
for k = 1:numel(ev)
    idx = ev(k)+sRel; if idx(1)<1 || idx(end)>nST, continue; end
    xx = double(m.Data.x(chans,idx))*lsb; xx = xx - median(xx(ri,:),1);
    prof = prof + filtfilt(bS,aS, xx(ci,:)')'; n = n+1;
end
prof = prof/n; sm = movmean(prof,5,1); v = mean(sm(:,cm),2);
v(1:5) = NaN; v(end-4:end) = NaN;
[~,iRad] = min(v); radCh = chOrder(iRad);
end


% ===================================================================
function [peakZ,swZ,fPeak,nCyc,ratioHi] = local_features(m,nST,fs,lsb,chOrder, ...
        pyrCh,radCh,refList,band,swBand,rip)
chans = unique([pyrCh radCh refList]); pidx = find(chans==pyrCh); sidx = find(chans==radCh); ri = ismember(chans,refList);
[bR,aR] = butter(3, band/(fs/2),'bandpass');
[bS,aS] = butter(3, swBand/(fs/2),'bandpass');
hi = [min(300,floor(0.45*fs)) min(400,floor(0.48*fs))]; [bH,aH] = butter(3, hi/(fs/2),'bandpass');
w = [-0.3 0.3]; sRel = round(w(1)*fs):round(w(2)*fs); t = sRel/fs;
em = t>=-0.05 & t<=0.05; bm = (t>=-0.25 & t<=-0.12) | (t>=0.12 & t<=0.25);
nfft = 2^nextpow2(round(0.08*fs));
N = size(rip,1);
peakZ = nan(N,1); swZ = peakZ; fPeak = peakZ; nCyc = peakZ; ratioHi = peakZ;
for k = 1:N
    idx = rip(k,2)+sRel; if idx(1)<1 || idx(end)>nST, continue; end
    xx = double(m.Data.x(chans,idx))*lsb; car = median(xx(ri,:),1);
    py = xx(pidx,:)-car; rd = xx(sidx,:)-car;
    fr = filtfilt(bR,aR,py); fh = filtfilt(bH,aH,py); fsw = filtfilt(bS,aS,rd);
    env = abs(hilbert(fr)); mu = median(env(bm)); sd = median(abs(env(bm)-mu))*1.4826;
    if sd==0 || isnan(sd), continue; end
    peakZ(k)  = max((env(em)-mu)/sd);
    ratioHi(k)= mean(fr(em).^2)/(mean(fh(em).^2)+eps);
    st = (rip(k,1)-rip(k,2))/fs; sp = (rip(k,3)-rip(k,2))/fs; dm = t>=st & t<=sp;
    y = fr(dm); nCyc(k) = sum(abs(diff(sign(y)))>0)/2;
    [pe,fg] = pwelch(py(em),[],[],nfft,fs); bb = fg>=120 & fg<=240; fg = fg(bb); pe = pe(bb);
    [~,ip] = max(pe); fPeak(k) = fg(ip);
    se = abs(hilbert(fsw)); mus = median(se(bm)); sds = median(abs(se(bm)-mus))*1.4826;
    if sds>0, swZ(k) = max((se(em)-mus)/sds); end
end
end


% ===================================================================
function idx = local_epochidx(t, ep)
% Return, for each time t, the index of the epoch [start stop] it falls in (0 if none).
idx = zeros(numel(t),1);
for e = 1:size(ep,1)
    idx(t>=ep(e,1) & t<=ep(e,2)) = e;
end
end


% ===================================================================
function local_qcfig(m,nST,fs,lsb,chOrder,iPyr,iRad,pyrCh,radCh,band,swBand,ev,R,outpng)
nC = numel(chOrder);
[bR,aR] = butter(3, band/(fs/2),'bandpass'); [bS,aS] = butter(3, swBand/(fs/2),'bandpass');
w = [-0.1 0.1]; sRel = round(w(1)*fs):round(w(2)*fs); t = sRel/fs;
chans = chOrder; ri = 1:nC;
sumDepth = zeros(nC,numel(sRel)); mPyrEnv = zeros(1,numel(sRel)); mSW = mPyrEnv; n = 0;
evs = ev; if numel(evs)>800, evs = evs(round(linspace(1,numel(evs),800))); end
for k = 1:numel(evs)
    idx = evs(k)+sRel; if idx(1)<1 || idx(end)>nST, continue; end
    xx = double(m.Data.x(chans,idx))*lsb; xx = xx - median(xx,1);
    sumDepth = sumDepth + xx;
    mPyrEnv = mPyrEnv + abs(hilbert(filtfilt(bR,aR,xx(iPyr,:))));
    mSW = mSW + filtfilt(bS,aS,xx(iRad,:));
    n = n+1;
end
sumDepth = sumDepth/n; mPyrEnv = mPyrEnv/n; mSW = mSW/n;

h = figure('Name',[R.sessionName ' QC'],'Position',[100 100 1100 400],'Visible','off');
subplot(1,3,1)
plot(t*1e3,mPyrEnv,'r'); hold on; plot(t*1e3,mSW,'b')
legend('pyr ripple env','rad SW'); xlabel('time from peak (ms)'); ylabel('\muV')
title(sprintf('%s  pyr %d / rad %d', R.sessionName, pyrCh, radCh),'Interpreter','none')
subplot(1,3,2)
imagesc(t,1:nC,sumDepth); axis xy; hold on; plot([0 0],[1 nC],'w:')
xlabel('time (s)'); ylabel('channel (depth)'); title('mean LFP across depth'); colormap('parula'); colorbar
subplot(1,3,3)
rt = R.peakTime; durTot = nST/fs; edges = 0:30:durTot;
stairs(edges(1:end-1)/60, histcounts(rt,edges)/30, 'k'); xlim([0 durTot/60])
xlabel('time (min)'); ylabel('events/s')
title(sprintf('band [%d %d] Hz | n=%d', band(1), band(2), R.nCurated))
try, saveas(h, outpng); catch, end
close(h)
end
