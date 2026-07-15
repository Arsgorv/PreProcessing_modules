function z = zscore_sliding(x, winSize)
%ZSCORE_SLIDING  Sliding-window z-score with NaN handling.
%
%   z = zscore_sliding(x, winSize)
%
%   x       : vector or matrix (NaNs allowed)
%   winSize : window length in samples (odd number recommended)
%
%   z       : z-scored data, same size as x
%
% Example:
%   z = zscore_sliding(signal, 100);

    if nargin < 2
        error('Need input data and window size');
    end
    
    % Mean & std over moving window, omitting NaNs
    mu  = movmean(x, winSize, 'omitnan');
    sig = movstd(x, winSize, 'omitnan');
    
    % Avoid divide-by-zero
    sig(sig == 0) = NaN;
    
    % Compute z-score
    z = (x - mu) ./ sig;
end
