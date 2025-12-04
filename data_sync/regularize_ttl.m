function [t_reg_s, info] = regularize_ttl(t_raw_s, Nexpected, name, plt)
% regularize_ttl  Insert missing TTLs and remove false ones
%
% INPUT
%   t_raw_s   : column vector of TTL times in seconds (absolute or relative)
%   Nexpected : expected number of events (eg N frames)
%   name      : string for messages ('fUS','video','baphy',...)
%   plt       : 0/1, make quick debug figure
%
% OUTPUT
%   t_reg_s   : regularized TTL times in seconds
%   info      : struct with debug info

if nargin < 3 || isempty(name)
    name = 'ttl';
end
if nargin < 4
    plt = 0;
end

t = t_raw_s(:);

if numel(t) < 2
    error('regularize_ttl:%s:NotEnough', name);
end

% -------------------------------------------------------------------------
% 1) Estimate the "frame ITI" from short intervals only (ignore huge gaps)
% -------------------------------------------------------------------------
iti0 = diff(t);
iti0 = iti0(:);

% sort and take lower 60–70% as "short" (frame-to-frame) intervals
iti_sorted = sort(iti0);
if numel(iti_sorted) >= 10
    n_short = round(0.7 * numel(iti_sorted));
    short_iti = iti_sorted(1:n_short);
else
    short_iti = iti_sorted;
end

medITI = median(short_iti);

% thresholds relative to frame ITI
big_thr        = 1.3 * medITI;   % gap > big_thr & < huge_thr => missing event(s)
small_thr      = 0.6 * medITI;   % gap < small_thr            => false event
huge_gap_thr   = 5.0 * medITI;   % gap >= huge_gap_thr        => "true" segment boundary

max_iter = 30;
iter = 0;

while numel(t) ~= Nexpected && iter < max_iter
    iter = iter + 1;
    
    iti = diff(t);
    
    % 1) insert missing events in moderately big gaps (but not across true
    %    segment boundaries)
    big_idx = find(iti > big_thr & iti < huge_gap_thr);
    if ~isempty(big_idx)
        extra = 0;
        for k = 1:numel(big_idx)
            idx = big_idx(k) + extra;
            new_t = (t(idx) + t(idx+1)) / 2;
            t = [t(1:idx); new_t; t(idx+1:end)];
            extra = extra + 1;
        end
    end
    
    % 2) remove false events (too short gaps)
    iti = diff(t);
    small_idx = find(iti < small_thr);
    if ~isempty(small_idx)
        t(small_idx + 1) = [];
    end
    
    % 3) gentle trim of the end if we overshoot a tiny bit
    diffN = numel(t) - Nexpected;
    if diffN > 0 && diffN <= 5
        t(end-diffN+1:end) = [];
    elseif diffN > 5
        % something is off; don't massacre data, just break
        warning('regularize_ttl:%s:TooManyAfterReg (diff=%d)', name, diffN);
        break
    end
end

t_reg_s = t;

info = struct;
info.name          = name;
info.N_raw         = numel(t_raw_s);
info.N_expected    = Nexpected;
info.N_final       = numel(t);
info.n_iter        = iter;
info.medITI_base   = medITI;
info.big_thr       = big_thr;
info.small_thr     = small_thr;
info.huge_gap_thr  = huge_gap_thr;

if numel(t) > 1
    info.medITI_final = median(diff(t));
else
    info.medITI_final = NaN;
end

if info.N_final ~= Nexpected
    warning('regularize_ttl:%s:Mismatch final=%d expected=%d', ...
        name, info.N_final, info.N_expected);
end

if plt
    figure;
    subplot(2,1,1);
    plot(t_raw_s, zeros(size(t_raw_s)),'k.');
    hold on;
    plot(t_reg_s, 0.1*ones(size(t_reg_s)),'r.');
    ylim([-0.2 0.2])
    title(sprintf('%s TTL regularization: raw=%d, final=%d, expected=%d', ...
        name, info.N_raw, info.N_final, info.N_expected));
    legend({'raw','reg'});
    
    subplot(2,1,2);
    histogram(diff(t_raw_s),50);
    hold on;
    histogram(diff(t_reg_s),50);
    xlabel('ITI (s)');
    legend({'raw','reg'});
end
end
