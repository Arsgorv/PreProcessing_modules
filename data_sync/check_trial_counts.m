function check_trial_counts(n_oe, n_baphy, n_csv, datapath, phase)
% check_trial_counts
% Compare trial counts from different sources.
%
% INPUT
%   n_oe    : #trials inferred from OE Baphy TTLs (can be empty)
%   n_baphy : #trials in Baphy metadata (B.n_trials)
%   n_csv   : #trials inferred from csv Baphy TTLs (can be empty)
%   datapath: session path
%   phase   : string label ('RP_passive', 'post-sync', etc.)

if nargin < 5
    phase = '';
end
if isempty(datapath), datapath = 'unknown'; end

% Represent empties nicely
if isempty(n_oe),    s_oe    = '?'; else, s_oe    = num2str(n_oe);    end
if isempty(n_baphy), s_baphy = '?'; else, s_baphy = num2str(n_baphy); end
if isempty(n_csv),   s_csv   = '?'; else, s_csv   = num2str(n_csv);   end

msg = sprintf('[%s] trials in %s: OE=%s, Baphy=%s, CSV=%s', ...
    phase, datapath, s_oe, s_baphy, s_csv);
disp(msg);

% Compare OE vs CSV
if ~isempty(n_oe) && ~isempty(n_csv) && n_csv ~= n_oe
    warning('[%s] Trial count mismatch OE vs CSV in %s: OE=%d, CSV=%d', ...
        phase, datapath, n_oe, n_csv);
end

% Compare OE vs Baphy
if ~isempty(n_oe) && ~isempty(n_baphy) && n_baphy ~= n_oe
    warning('[%s] Trial count mismatch OE vs Baphy in %s: OE=%d, Baphy=%d', ...
        phase, datapath, n_oe, n_baphy);
end

% Compare Baphy vs CSV
if ~isempty(n_baphy) && ~isempty(n_csv) && n_baphy ~= n_csv
    warning('[%s] Trial count mismatch Baphy vs CSV in %s: Baphy=%d, CSV=%d', ...
        phase, datapath, n_baphy, n_csv);
end

end
