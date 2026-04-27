function cfg = np_lfp_config(datapath)
% Per-animal NP-LFP setup. Add new animals here; return [] for animals
% without NP-in-HPC so the block is a no-op.
cfg = [];
if contains(datapath, 'Tvorozhok')
    cfg.probe         = 'A';
    cfg.channels      = 350;        % or [deep mid1 mid2 sup]
    cfg.theta_channel = 350;
elseif contains(datapath, 'Mochi')
    cfg.probe         = 'A';
    cfg.channels      = 300;
    cfg.theta_channel = 300;
end
end