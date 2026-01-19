function trig_csv = regularize_all_ttl_csv(trig_csv, Baphy, datapath)
% regularize_all_ttl_csv
% For csv-only sessions (e.g. Chabichou):
%   - Regularise fUS TTL train (fus_time_s) strongly (grid-based)
%   - Keep Baphy TTLs structurally intact (no cross-trial interpolation),
%     but run generic clean-up inside bursts if needed.

if ~isstruct(trig_csv)
    return
end

info = struct;

%% 1) fUS TTLs (csv, only Exp)
if isfield(trig_csv,'fus_time_s') && ~isempty(trig_csv.fus_time_s)
    [t_reg, info_fus] = regularize_ttl_sequence(trig_csv.fus_time_s(:), 'csv-fUS', datapath);
    trig_csv.fus_time_s = t_reg;
    info.fus = info_fus;

    % check against Exp frames in exp_info (if available)
    try
        fus_dir = fullfile(datapath,'fUS');
        exp_file = fullfile(fus_dir,'exp_info.mat');
        if exist(exp_file,'file')
            S = load(exp_file);
            exp_info = S.exp_info;

            if isfield(exp_info,'PreExp')
                n_exp = exp_info.Exp.size(3);
            elseif isfield(exp_info,'size')
                n_exp = exp_info.size{2}(3);
            else
                n_exp = [];
            end

            if ~isempty(n_exp)
                n_fus = numel(trig_csv.fus_time_s);
                if n_fus ~= n_exp
                    warning('[csv-fUS] Exp frames mismatch in %s: TTL=%d, fUS=%d', ...
                        datapath, n_fus, n_exp);
                end
            end
        end
    catch ME
        warning('[csv-fUS] exp_info Exp frame check failed in %s: %s', datapath, ME.message);
    end
else
    info.fus = [];
end

%% 2) Baphy TTLs – as before (no cross-trial interpolation)
if isfield(trig_csv,'baphy_time_s') && ~isempty(trig_csv.baphy_time_s)
    [t_reg_b, info_baphy] = regularize_ttl_sequence(trig_csv.baphy_time_s(:), 'csv-Baphy', datapath);
    trig_csv.baphy_time_s = t_reg_b;

    t_baphy = trig_csv.baphy_time_s(:);
    try
        [~, ~, first_idx, last_idx] = detect_trial_onsets_from_baphy_ttl(datapath, t_baphy, B.n_trials);
    catch
        [~, ~, first_idx, last_idx] = detect_trial_onsets_from_baphy_ttl(datapath, t_baphy, []);
    end

    if ~isempty(first_idx)
        n_per_trial = last_idx - first_idx + 1;
        info_baphy.n_per_trial = n_per_trial;
        info_baphy.n_mode      = mode(n_per_trial);
    end

    info.baphy = info_baphy;
else
    info.baphy = [];
end

%{
    if contains(datapath, ['Chabichou' filesep '20210406_1_m_F'])
        b_trigs = [b_trigs(1:162)' b_trigs(162)+1000 b_trigs(162)+2000 b_trigs(163:end)']';
    elseif contains(datapath, ['Chabichou' filesep '20210413_1_m_I'])
        b_trigs(31) = [];
        b_trigs(31) = [];
        b_trigs(31) = [];
    end
%}

trig_csv.info_csv_regularize = info;

end
