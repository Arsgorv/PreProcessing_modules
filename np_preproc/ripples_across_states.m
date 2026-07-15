%% ripples_across_states.m
%  Quantify ripple specificity across vigilance states using intervalSets.
%  Requires TStoolbox on path (intervalSet, ts, Restrict, Start, End).
%
%  Provide in the workspace before running:
%    R          : ripple result from detect_ripples (peakTime in s, startStop
%                 in s, duration, features). Or load it (see below).
%    states     : cell array of intervalSet, one per vigilance state
%    stateNames : cell array of char, names matching 'states'
%    u          : tsd/intervalSet time units (time_sec * u). MUST match your
%                 scoring intervalSets. TStoolbox default = 1e4.
%
%  Computes per state: ripple count, duration, rate (events/s), share of all
%  ripples, share of scored time, and enrichment (ripple% / time%). >1 =
%  ripples over-represented in that state.

%% --- config / inputs ---
u = 1e4;                              % <-- MATCH your intervalSet time units!

% Example wiring (uncomment & adapt to how you store scoring):
% S = load(fullfile(R.datapath,'SleepState.mat'));        % your scoring file
% states     = {S.Wake, S.NREM, S.REM};
% stateNames = {'Wake','NREM','REM'};

assert(exist('R','var')==1,        'Load a ripple result into R first.');
assert(exist('states','var')==1,   'Define cell array ''states'' (intervalSets).');
assert(exist('stateNames','var')==1,'Define cell array ''stateNames''.');

pk = R.peakTime(:);                  % ripple peak times (s)
rippleTS = ts(pk*u);                 % peaks as ts (for Restrict)
rippleIS = intervalSet(R.startStop(:,1)*u, R.startStop(:,2)*u);   % usable downstream

%% --- unit sanity check ---
maxEndSec = max(cellfun(@(is) max(End(is)), states)) / u;
if isfield(R,'recDurSec') && maxEndSec > 1.5*R.recDurSec
    warning(['State end times reach %.0f s but the recording is %.0f s. ', ...
        'Your ''u'' (tsd units) is probably wrong.'], maxEndSec, R.recDurSec);
end

%% --- per-state statistics ---
nS  = numel(states);
cnt = zeros(nS,1); durS = zeros(nS,1);
stateOfRipple = zeros(numel(pk),1);          % per-ripple state index (0 = none)
for s = 1:nS
    seg = [Start(states{s}) End(states{s})] / u;     % state intervals (s)
    durS(s) = sum(seg(:,2)-seg(:,1));
    inAny = false(numel(pk),1);
    for r = 1:size(seg,1)
        inAny = inAny | (pk>=seg(r,1) & pk<=seg(r,2));
    end
    cnt(s) = sum(inAny);
    stateOfRipple(inAny) = s;
end
rate   = cnt ./ max(durS,eps);               % events / s
propR  = cnt / numel(pk);                     % share of all ripples
propT  = durS / sum(durS);                    % share of scored time
enrich = propR ./ max(propT,eps);             % >1 = over-represented

%% --- report ---
fprintf('\nRipples across states (n=%d total, %d outside all states):\n', ...
    numel(pk), numel(pk)-sum(cnt));
fprintf('%-10s %8s %10s %9s %10s %9s\n','state','n','dur(s)','rate/s','%%ripples','enrich');
for s = 1:nS
    fprintf('%-10s %8d %10.1f %9.3f %9.1f%% %9.2f\n', ...
        stateNames{s}, cnt(s), durS(s), rate(s), 100*propR(s), enrich(s));
end

%% --- figures ---
figure('Name','ripples across states')
subplot(1,3,1)
bar(rate); set(gca,'xticklabel',stateNames); ylabel('ripple rate (events/s)'); title('rate by state')
subplot(1,3,2)
bar([propT propR]); set(gca,'xticklabel',stateNames); ylabel('proportion')
legend('time','ripples'); title('time vs ripple share')
subplot(1,3,3)
bar(enrich); hold on; plot(xlim,[1 1],'k:'); set(gca,'xticklabel',stateNames)
ylabel('enrichment (ripple%% / time%%)'); title('state specificity')

% feature distributions by state (common bins per feature)
if isfield(R,'features')
    feats = {R.features.peakZ(:),'ripple z'; R.duration(:)*1e3,'duration (ms)'; ...
             R.features.peakFreq(:),'peak freq (Hz)'};
    figure('Name','ripple features by state')
    for fi = 1:size(feats,1)
        v = feats{fi,1}; edges = linspace(min(v),max(v),21); ctr = edges(1:end-1)+diff(edges)/2;
        subplot(1,3,fi); hold on
        for s = 1:nS
            vs = v(stateOfRipple==s); if isempty(vs), continue; end
            n = histcounts(vs,edges); plot(ctr, n/sum(n), 'LineWidth',1.2)
        end
        xlabel(feats{fi,2}); ylabel('norm. count'); if fi==1, legend(stateNames); end
    end
end

%% --- per-state ripple sets for downstream (TStoolbox) ---
% e.g. ripples during NREM, and any signal restricted to them:
%   nremRippleTS = Restrict(rippleTS, states{strcmp(stateNames,'NREM')});
%   lfpDuringNREMripples = Restrict(myLFPtsd, intersect(rippleIS, states{...}));
rippleByState = cell(nS,1);
for s = 1:nS
    rippleByState{s} = Restrict(rippleTS, states{s});
end
