function check_trial_counts(n_oe, n_baphy, n_csv, datapath, phase)

msg = sprintf('[%s] trials: OE=%d, Baphy=%d, CSV=%d', phase, n_oe, n_baphy, n_csv);
disp(msg);

if ~isempty(n_csv) && n_csv ~= n_oe
    warning('Trial count mismatch OE vs CSV in %s: OE=%d, CSV=%d', datapath, n_oe, n_csv);
end

if ~isempty(n_baphy) && n_baphy ~= n_oe
    warning('Trial count mismatch OE vs Baphy in %s: OE=%d, Baphy=%d', datapath, n_oe, n_baphy);
end
end