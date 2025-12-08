function trigOE = regularize_all_ttl(trigOE, Baphy, datapath)
% regularize_all_ttl
% Regularise OE TTLs (fUS / video / Baphy, etc.) using regularize_ttl_sequence.
%
% For Baphy we do NOT interpolate across trial pauses; only stats.
% Trial-wise interpolation (if ever needed) should be done explicitly.

if nargin < 3
    datapath = '';
end

channel_list = {'fus','baphy','face_cam','eye_cam'};

for c = 1:numel(channel_list)
    fname = channel_list{c};
    if isfield(trigOE, fname) && isfield(trigOE.(fname),'t_raw_s') ...
            && ~isempty(trigOE.(fname).t_raw_s)

        t_raw = trigOE.(fname).t_raw_s(:);
        label = ['OE-' fname];

        [t_reg, info] = regularize_ttl_sequence(t_raw, label, datapath);

        % store regularised version in seconds
        trigOE.(fname).t_s = t_reg;
        trigOE.(fname).regularize_info = info;

        % For channels where we really want to use the regularised times as
        % "raw", update t_raw_s too (especially fUS / video)
        if contains(lower(fname),'fus') || contains(lower(fname),'video')
            trigOE.(fname).t_raw_s = t_reg;
        end
    end
end

%% fUS frame-count sanity check vs exp_info (OE case only)
if isfield(trigOE,'fus') && isfield(trigOE.fus,'t_s') && ~isempty(trigOE.fus.t_s)
    try
        fus_dir = fullfile(datapath,'fUS');
        exp_file = fullfile(fus_dir,'exp_info.mat');
        if exist(exp_file,'file')
            S = load(exp_file);
            exp_info = S.exp_info;

            if isfield(exp_info,'PreExp')
                n_pre_exp  = exp_info.PreExp.size(3);
                n_exp      = exp_info.Exp.size(3);
                n_post_exp = exp_info.PostExp.size(3);
            elseif isfield(exp_info,'size')
                n_pre_exp  = exp_info.size{1}(3);
                n_exp      = exp_info.size{2}(3);
                n_post_exp = exp_info.size{3}(3);
            else
                n_pre_exp = [];
                n_exp     = [];
                n_post_exp= [];
            end

            t = trigOE.fus.t_s(:);
            dt = diff(t);
            gap_idx = find(dt > 3); % same phase-gap threshold
            seg_start = [1; gap_idx + 1];
            seg_end   = [gap_idx; numel(t)];
            n_seg = numel(seg_start);
            n_fus_seg = seg_end - seg_start + 1;

            if n_seg >= 3 && ~isempty(n_pre_exp) && ~isempty(n_exp) && ~isempty(n_post_exp)
                n_pre_fus  = n_fus_seg(1);
                n_exp_fus  = n_fus_seg(2);
                n_post_fus = n_fus_seg(3);

                if n_pre_fus ~= n_pre_exp
                    warning('[OE-fUS] PreExp frames mismatch in %s: TTL=%d, fUS=%d', ...
                        datapath, n_pre_fus, n_pre_exp);
                end
                if n_exp_fus ~= n_exp
                    warning('[OE-fUS] Exp frames mismatch in %s: TTL=%d, fUS=%d', ...
                        datapath, n_exp_fus, n_exp);
                end
                if n_post_fus ~= n_post_exp
                    warning('[OE-fUS] PostExp frames mismatch in %s: TTL=%d, fUS=%d', ...
                        datapath, n_post_fus, n_post_exp);
                end
            else
                warning('[OE-fUS] Could not clearly segment 3 fUS blocks in %s (n_seg=%d).', ...
                    datapath, n_seg);
            end
        end
    catch ME
        warning('[OE-fUS] exp_info frame check failed in %s: %s', datapath, ME.message);
    end
end

% For Baphy specifically, keep t_raw_s unchanged except for potential
% future use; Edel trial logic uses trigOE.baphy.t_raw_s (not interpolated).

end
