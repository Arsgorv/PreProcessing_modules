function LH = dlc_get_lh(D, bpIdx)
% Return likelihood for a single-point bodypart (or first if multiple).
if isempty(bpIdx)
    LH = nan(size(D,1),1);
    return
end
bp = bpIdx(1);
col = 2 + (bp-1)*3 + 3;
LH = D(:, col);
LH = LH(:);
end