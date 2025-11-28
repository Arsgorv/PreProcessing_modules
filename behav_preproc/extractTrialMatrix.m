function [mat, keepIdx] = extractTrialMatrix( ...
            tVec, y, onsets, baselineN, trialN, smoothWin)
%EXTRACTTRIALMATRIX  Cut a continuous signal into baseline-corrected trials
%
%   [MAT, KEEPIDX] = extractTrialMatrix(tVec, y, onsets, ...
%                                       baselineN, trialN, smoothWin)
%
% Inputs
%   tVec       : time stamps of y   (same units as ONSETS, typically seconds)
%   y          : signal vector      (length = numel(tVec))
%   onsets     : trial-start times  (Nx1 double)
%   baselineN  : # samples BEFORE onset for baseline mean subtraction
%   trialN     : # samples to keep  (includes the onset sample itself)
%   smoothWin  : smoothing window (samples); set = 0 to skip smoothing
%
% Outputs
%   MAT        : nTrials × trialN matrix (NaN for rejected trials)
%   KEEPIDX    : logical length(onsets) flag of successfully extracted trials
%
% Notes
%   - Uses runmean smoothing so outliers don’t distort the baseline.
%   - A trial is rejected if the full baseline or trial window would
%     fall outside the range of the signal.
%
% Arsenii Goriachenkov
% July-2025
% LSP - MOBS
% -------------------------------------------------------------------------

nTrials   = numel(onsets);
mat       =  NaN(nTrials, trialN);
keepIdx   = false(nTrials,1);

% ---- pre-compute smoothing parameter for runmean ------------------------
if smoothWin > 1
    % smoothWin is a *total* window length in samples (e.g. 15)
    % runmean expects half-window M so that window = 2*M+1
    M = floor((smoothWin-1)/2);
else
    M = 0;   % no smoothing
end

for k = 1:nTrials
    t0   = onsets(k);
    idx0 = find(tVec >= t0, 1, 'first');     % first sample at/after onset
    if isempty(idx0), continue; end

    idxBaseline = idx0-baselineN : idx0-1;
    idxTrial    = idx0 : idx0+trialN-1;

    % Reject if window exceeds signal boundaries
    if idxBaseline(1) < 1 || idxTrial(end) > numel(y), continue; end

    baseline   = mean(y(idxBaseline), 'omitnan');
    raw        = y(idxTrial) - baseline;     % ? from baseline
    
    % optional NaN handling – runmean cannot ignore NaNs
    if any(isnan(raw))
        raw = fillmissing(raw,'linear','EndValues','nearest');
    end
    
    % --- smooth with runmean --------------------------------------------
    if M > 0
        smoothed = runmean(raw, M);    % runs along the row vector
    else
        smoothed = raw;
    end

    mat(k,:)  = smoothed;
    keepIdx(k)= true;
end
end