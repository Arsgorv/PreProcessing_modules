function LFP = load_channel_as_tsd(outMatFile, ch)
% load_channel_as_tsd
% Load one channel from an LFP matrix file produced by Master_LFP_NP_preproc (matrix mode).
% Uses matfile to avoid loading the whole Y matrix.

if exist(outMatFile,'file') ~= 2
    error('load_channel_as_tsd:NoFile','Missing file: %s', outMatFile);
end

M = matfile(outMatFile);

if ~isprop(M,'chUse') || ~isprop(M,'Y')
    error('load_channel_as_tsd:BadFile','Missing chUse/Y in %s', outMatFile);
end

chUse = M.chUse;
idx = find(chUse == ch, 1, 'first');
if isempty(idx)
    error('load_channel_as_tsd:MissingChannel','Channel %d not present in %s', ch, outMatFile);
end

if isprop(M,'t0_ts') && isprop(M,'dt_ts') && isprop(M,'nOut')
    t0_ts = double(M.t0_ts);
    dt_ts = double(M.dt_ts);
    nOut  = double(M.nOut);
    t_ts = t0_ts + (0:nOut-1)' * dt_ts;
else
    error('load_channel_as_tsd:NoTime','Missing t0_ts/dt_ts/nOut in %s', outMatFile);
end

y = M.Y(:, idx);
LFP = tsd(t_ts, double(y));
end