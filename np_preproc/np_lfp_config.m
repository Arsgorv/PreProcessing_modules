function cfg = np_lfp_config(datapath)
% Per-animal NP-LFP setup. Add new animals here; return [] for animals
% without NP-in-HPC so the block is a no-op.
%
% Fields:
%   probe          : 'A' | 'B'   NP probe carrying the HPC LFP
%   channels       : 1-based channel(s) in the Probe*-LFP stream to export
%   theta_channel  : channel written to ChannelsToAnalyse/ThetaREM.mat
%   master_ttl_chan: LFP channel (in LFPData/) carrying the OneBox sync pulse on
%                    the Acquisition Board. Used for analog OneBox<->master sync.
%                    Must match cfg.OneBox in get_trigger_config for this animal.
%                    Leave empty to fall back to get_trigger_config (may prompt).
%   onebox_adc_idx : 1-based channel in the OneBox-ADC stream carrying the same
%                    pulse (default 3 = "ADC2").
cfg = [];
if contains(datapath, 'Tvorozhok')
    cfg.probe           = 'A';
    cfg.channels        = 350;        % or [deep mid1 mid2 sup]
    cfg.theta_channel   = 350;
    cfg.master_ttl_chan = 91;         % Tvorozhok experiment: cfg.OneBox = 91
    cfg.onebox_adc_idx  = 3;
elseif contains(datapath, 'Mochi')
    cfg.probe           = 'A';
    cfg.channels        = 125;
    cfg.theta_channel   = 125;
    cfg.master_ttl_chan = 56;         % Mochi experiment: cfg.OneBox = 56
    cfg.onebox_adc_idx  = 3;
end
end