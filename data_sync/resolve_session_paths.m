function P = resolve_session_paths(datapath)
% resolve_session_paths
% Canonical resolver for variable session layouts (Active + Passive).

P.datapath = datapath;

P.video_dir = fullfile(datapath, 'video');
P.stim_dir  = fullfile(datapath, 'stim');
P.baphy_dir = fullfile(datapath, 'baphy');
P.ephys_dir = fullfile(datapath, 'ephys');

P.lfp_dir = '';
if exist(fullfile(datapath, 'LFPData'), 'dir')
    P.lfp_dir = fullfile(datapath, 'LFPData');
elseif exist(fullfile(datapath, 'ephys', 'LFPData'), 'dir')
    P.lfp_dir = fullfile(datapath, 'ephys', 'LFPData');
end

P.has_video = exist(P.video_dir, 'dir') ~= 0;
P.has_stim  = exist(P.stim_dir,  'dir') ~= 0;
P.has_baphy = exist(P.baphy_dir, 'dir') ~= 0;
P.has_lfp   = ~isempty(P.lfp_dir) && exist(P.lfp_dir, 'dir') ~= 0;
end
