function [t_reg, info] = regularize_ttl_sequence(t_raw, label, datapath)

if nargin < 2 || isempty(label)
    label = 'TTL';
end
if nargin < 3
    datapath = '';
end

t_reg = [];
info = struct('n_raw',0,'n_final',0,'median_iti',NaN, ...
              'n_removed',0,'n_inserted',0,'removed_idx',[]);

if isempty(t_raw)
    return
end

t = sort(t_raw(:));
info.n_raw = numel(t);

if numel(t) < 2
    t_reg = t;
    info.n_final = numel(t);
    return
end

lbl = lower(label);
iti = diff(t);
medITI = median(iti);
info.median_iti = medITI;

%% === Case 1: fUS (OE or csv) ============================================
if contains(lbl,'fus')
    if medITI <= 0
        warning('[%s] Non-positive median ITI for fUS in %s', label, datapath);
        t_reg = t;
        info.n_final = numel(t);
        return
    end

    % 1) find big gaps between phases (PreExp / Exp / PostExp)
    phase_gap_thr = 3; % seconds
    gap_idx = find(diff(t) > phase_gap_thr);
    seg_start = [1; gap_idx + 1];
    seg_end   = [gap_idx; numel(t)];
    n_seg = numel(seg_start);

    n_removed_total  = 0;
    n_inserted_total = 0;
    t_all = [];

    % 2) regularise each segment separately
    for s = 1:n_seg
        idx_range = seg_start(s):seg_end(s);
        t_seg = t(idx_range);

        if numel(t_seg) < 2
            t_all = [t_all; t_seg]; %#ok<AGROW>
            continue
        end

        iti_seg = diff(t_seg);
        medITI_seg = median(iti_seg);

        if medITI_seg <= 0
            t_all = [t_all; t_seg]; %#ok<AGROW>
            continue
        end

        lower_thr = 0.7 * medITI_seg;  % too close -> extra
        upper_thr = 1.4 * medITI_seg;  % too far  -> missing

        n_removed_seg  = 0;
        n_inserted_seg = 0;

        max_iter = 10;
        for iter = 1:max_iter
            changed = false;

            % --- remove extra close pulses ---
            while true
                iti_seg = diff(t_seg);
                idx_close = find(iti_seg < lower_thr, 1, 'first');
                if isempty(idx_close)
                    break
                end
                t_seg(idx_close+1) = [];
                n_removed_seg = n_removed_seg + 1;
                changed = true;
            end

            % --- insert missing pulses in long gaps ---
            iti_seg = diff(t_seg);
            idx_long = find(iti_seg > upper_thr);
            if ~isempty(idx_long)
                new_all = [];
                for k = 1:numel(idx_long)
                    i = idx_long(k);
                    dt_gap = iti_seg(i);
                    n_missing = round(dt_gap / medITI_seg) - 1;
                    if n_missing <= 0
                        continue
                    end
                    new_times = t_seg(i) + (1:n_missing) * (dt_gap / (n_missing + 1));
                    new_all = [new_all; new_times(:)]; %#ok<AGROW>
                    n_inserted_seg = n_inserted_seg + n_missing;
                    changed = true;
                end
                if ~isempty(new_all)
                    t_seg = sort([t_seg; new_all]);
                end
            end

            if ~changed
                break
            end
        end

        n_removed_total  = n_removed_total + n_removed_seg;
        n_inserted_total = n_inserted_total + n_inserted_seg;

        t_all = [t_all; t_seg]; 
    end
    
    % 20251204 AG: i'm tired of finding general solution for random shit, so I'll apply exceptions here, sorry.
    if contains(datapath, ['Edel' filesep '20220419_1_m_C'])
        t_all(1) = [];
        n_removed_total = n_removed_total + 1;
    end

    info.n_removed  = n_removed_total;
    info.n_inserted = n_inserted_total;
    info.n_final    = numel(t_all);

    if n_removed_total > 0 || n_inserted_total > 0
        warning('[%s] fUS TTLs regularised in %s: removed %d, inserted %d (raw=%d, final=%d)', ...
            label, datapath, n_removed_total, n_inserted_total, info.n_raw, info.n_final);
    end

    t_reg = t_all;
    
    return
end

%% === Case 2: Baphy (no interpolation across trial pauses) ===============
if contains(lbl,'baphy')
    if medITI <= 0
        warning('[%s] Non-positive median ITI for Baphy in %s', label, datapath);
    end
    t_reg = t;
    info.n_final = numel(t);
    return
end

%% === Case 3: generic (video, cam, etc.) =================================
if medITI <= 0
    warning('[%s] Non-positive median ITI in %s', label, datapath);
    t_reg = t;
    info.n_final = numel(t);
    return
end

lower_thr = 0.5 * medITI;
upper_thr = 1.5 * medITI;

removed_idx = [];
changed = true;
while changed
    changed = false;
    iti = diff(t);
    bad = find(iti < lower_thr, 1, 'first');
    if ~isempty(bad)
        t(bad+1) = [];
        removed_idx = [removed_idx; bad+1]; %#ok<AGROW>
        changed = true;
    end
end
info.removed_idx = removed_idx;
info.n_removed   = numel(removed_idx);

iti = diff(t);
big_idx = find(iti > upper_thr);
to_add = [];

for k = 1:numel(big_idx)
    i = big_idx(k);
    dt_gap = iti(i);
    n_missing = round(dt_gap / medITI) - 1;
    if n_missing <= 0
        continue
    end
    new_times = t(i) + (1:n_missing) * (dt_gap / (n_missing + 1));
    to_add = [to_add; new_times(:)]; %#ok<AGROW>
end

if ~isempty(to_add)
    t = sort([t; to_add]);
end

info.n_inserted = numel(to_add);
info.n_final    = numel(t);

if info.n_removed > 0 || info.n_inserted > 0
    warning('[%s] TTL regularisation in %s: removed %d, inserted %d (raw=%d, final=%d)', ...
        label, datapath, info.n_removed, info.n_inserted, info.n_raw, info.n_final);
end

t_reg = t;

end
