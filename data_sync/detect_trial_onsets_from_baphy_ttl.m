function [t_trial_onset, t_trial_offset, first_idx, last_idx] = ...
    detect_trial_onsets_from_baphy_ttl(datapath, t_baphy, n_trials_expected)
% Detect trial ONSET and OFFSET times from Baphy TTL times.
% Assumes within-trial ITIs ~1 s, between-trial ITIs ~2 s or more.

t_trial_onset  = [];
t_trial_offset = [];
first_idx = [];
last_idx  = [];

if isempty(t_baphy)
    return
end

t_baphy = t_baphy(:);

if contains(datapath, 'Edel\20220419_1_m_C') && numel(t_baphy) > 1100
    t_baphy(1) = [];
end


n_ttl = numel(t_baphy);

if n_ttl == 1
    t_trial_onset  = t_baphy;
    t_trial_offset = t_baphy;
    first_idx = 1;
    last_idx  = 1;
    return
end

iti = diff(t_baphy);

% --- Robust short vs long ITI separation ---
if numel(iti) >= 5
    q_low  = quantile(iti, 0.7);   % end of "short" ITIs
    q_high = quantile(iti, 0.95);  % clear "long" ITIs

    short_iti = iti(iti <= q_low);
    long_iti  = iti(iti >= q_high);

    if ~isempty(short_iti) && ~isempty(long_iti)
        mu_short = mean(short_iti);
        mu_long  = mean(long_iti);
        thr = (mu_short + mu_long) / 2;
    else
        thr = 1.5 * median(iti);
    end
else
    thr = 1.5 * median(iti);
end

% Gap indices separate trials
gap_idx = find(iti > thr);

first_idx = [1; gap_idx + 1];
last_idx  = [gap_idx; n_ttl];

t_trial_onset  = t_baphy(first_idx);
t_trial_offset = t_baphy(last_idx);

% Optional trimming to expected number of trials
if nargin >= 2 && ~isempty(n_trials_expected) && n_trials_expected > 0
    n_detected = numel(t_trial_onset);
    if n_detected ~= n_trials_expected
        warning(['[detect_trial_onsets_from_baphy_ttl] Trials: expected %d, ' ...
                 'detected %d. Truncating to min.'], n_trials_expected, n_detected);
        n_keep = min(n_trials_expected, n_detected);
        first_idx     = first_idx(1:n_keep);
        last_idx      = last_idx(1:n_keep);
        t_trial_onset = t_trial_onset(1:n_keep);
        t_trial_offset= t_trial_offset(1:n_keep);
    end
end

end
