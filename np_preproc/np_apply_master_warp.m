function t_master_s = np_apply_master_warp(t_np_s, datapath, segName)
% Apply the saved master-clock warp for a segment to arbitrary NP-clock times
% (e.g. spike times exported by Kilosort).
%
% Inputs:
%   t_np_s   : vector of times in NP-stream clock (seconds since segment start)
%   datapath : session root
%   segName  : segment folder basename (must match what NP LFP pipeline used)
%
% Output:
%   t_master_s : same shape as t_np_s, in master-clock seconds since segment start.
%
% The warp is loaded from <datapath>/analysis/master_clock_sync_<seg>.mat.

mf = fullfile(datapath,'analysis',['master_clock_sync_' segName '.mat']);
if exist(mf,'file') ~= 2
    error('np_apply_master_warp:NoSync','No saved sync at %s', mf);
end
S = load(mf);

switch S.method
    case {'ttl','oe-synchronizer'}
        t_master_s = interp1(S.np_anchor_s, S.master_anchor_s, t_np_s, 'linear', 'extrap');
    case 'linear-ratio'
        t_master_s = t_np_s * S.warp_ratio;
    case {'identity-no-master','identity'}
        t_master_s = t_np_s;
    otherwise
        warning('Unknown sync method "%s"; identity', S.method);
        t_master_s = t_np_s;
end
end