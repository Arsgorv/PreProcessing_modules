function XY = dlc_get_xy(D, bpIdx, which)
% Return nFrames x nPts matrix of X (default) or Y for given bodypart indices.
% D: synchronized matrix (time, frame, triplets...)
% bpIdx: vector of bodypart indices in 1..nBp
% which: 'x' or 'y'
if nargin < 3
    which = 'x';
end
if isempty(bpIdx)
    XY = nan(size(D,1),0);
    return
end

switch lower(which)
    case 'x'
        offset = 1;
    case 'y'
        offset = 2;
    otherwise
        error('dlc_get_xy: which must be x or y');
end

cols = 2 + (bpIdx(:)-1)*3 + offset;
XY = D(:, cols);

% Ensure column vector for single-point parts
if size(XY,2) == 1
    XY = XY(:);
end
end
