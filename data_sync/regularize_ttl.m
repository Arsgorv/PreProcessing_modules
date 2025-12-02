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
    error('RP_regularize_ttl:%s:NotEnough', name);
end

max_iter = 30;
iter = 0;

while numel(t) ~= Nexpected && iter < max_iter
    iter = iter + 1;
    
    iti = diff(t);
    medITI = median(iti);
    
    big_thr   = 1.2 * medITI;   % gap > big_thr => missing event
    small_thr = 0.7 * medITI;   % gap < small_thr => false event
    
    % 1) insert missing events in big gaps
    big_idx = find(iti > big_thr);
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
    
    % 3) if still too many, trim end
    if numel(t) > Nexpected
        t = t(1:Nexpected);
    end
end

t_reg_s = t;

info = struct;
info.name       = name;
info.N_raw      = numel(t_raw_s);
info.N_expected = Nexpected;
info.N_final    = numel(t);
info.n_iter     = iter;
if numel(t) > 1
    info.medITI_final = median(diff(t));
else
    info.medITI_final = NaN;
end

if info.N_final ~= Nexpected
    warning('regularize_ttl:%s:Mismatch', name);
end

if plt
    figure;
    subplot(2,1,1);
    plot(t_raw_s, zeros(size(t_raw_s)),'k.');
    hold on;
    plot(t_reg_s, 0.1*ones(size(t_reg_s)),'r.');
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
