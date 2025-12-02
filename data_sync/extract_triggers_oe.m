function trigOE = extract_triggers_oe(datapath, cfg)
% extract_triggers_oe
% Detect TTL pulses for fUS, Baphy, cameras, etc from OpenEphys LFPData.
% Optionally regularize fUS and video TTLs using RP_regularize_ttl.
%
% INPUT
%   datapath : session folder containing LFPData
%   cfg      : struct with fields (set [] or NaN if not used)
%       fus_ch       : LFP channel index for fUS TTLs
%       baphy_ch     : LFP channel index for Baphy TTLs
%       face_ch      : LFP channel index for face-cam TTLs
%       eye_ch       : LFP channel index for eye-cam TTLs
%       respi        : LFP channel index for respiration TTLs (Kosichka only)
%       heart        : LFP channel index for heart TTLs (Kosichka only)
%       OneBox       : LFP channel index for NP OneBox or similar, if any
%
%       N_fus_frames : expected # fUS frames (optional, for regularisation)
%       N_face_frames: expected # face frames (optional)
%       N_eye_frames : expected # eye frames (optional)
%
%       plt          : 0/1 – make quick sanity plots (default 0)
%
% OUTPUT
%   trigOE : struct with subfields:
%     .Fs_lfp           : sampling rate of LFP (Hz, estimated)
%     .time_lfp_ts      : time vector in ts units (1e-4 s)
%     .time_lfp_s       : time vector in seconds
%
%     For each TTL channel present (fus, baphy, face, eye, respi, heart, onebox):
%        .<name>.ch        : channel index
%        .<name>.peak_idx  : sample indices of TTL peaks
%        .<name>.peak_val  : peak values
%        .<name>.t_raw_ts  : times of peaks (ts units)
%        .<name>.t_raw_s   : times of peaks (seconds)
%        .<name>.t_reg_s   : regularized times (if N_expected provided)
%        .<name>.reg_info  : struct from RP_regularize_ttl
%
% NOTE
%   This function assumes that each LFP file LFP<ch>.mat contains a tsd
%   variable named LFP.

if ~isfield(cfg,'plt'), cfg.plt = 0; end

trigOE = struct;

% figure out animal name just for logging
animals = {'Edel','Chabichou','Kosichka','Tvorozhok','Shropshire','Brynza','Labneh'};
Session_params.animal = 'unknown';
for i = 1:numel(animals)
    if contains(datapath, animals{i})
        Session_params.animal = animals{i};
        break
    end
end

[~, sess_name] = fileparts(datapath);
disp(['[extract_triggers_oe] Animal: ' Session_params.animal ', session: ' sess_name]);

% helper to detect TTL peaks given a channel
ttl_names = {'fus','baphy','face','eye','respi','heart','onebox'};
ch_fields = {'fus_ch','baphy_ch','face_ch','eye_ch','respi','heart','OneBox'};

Fs_lfp = [];
time_lfp_ts = [];
time_lfp_s  = [];

for k = 1:numel(ttl_names)
    ttl_name = ttl_names{k};
    ch_field = ch_fields{k};
    if ~isfield(cfg, ch_field), continue; end
    ch = cfg.(ch_field);
    
    if isempty(ch) || isnan(ch)
        continue
    end
    
    lfp_file = fullfile(datapath,'ephys','LFPData', ['LFP' num2str(ch) '.mat']);
    if ~exist(lfp_file, 'file')
        warning('[extract_triggers_oe] LFP file not found for %s channel %d: %s', ttl_name, ch, lfp_file);
        continue
    end
    
    S = load(lfp_file);  % expects variable LFP
    if ~isfield(S, 'LFP')
        warning('[extract_triggers_oe] No LFP variable in %s', lfp_file);
        continue
    end
    LFP = S.LFP;
    
    t_ts = Range(LFP);      % ts units (1e-4 s)
    v    = Data(LFP);
    
    if isempty(time_lfp_ts)
        time_lfp_ts = t_ts(:);
        time_lfp_s  = time_lfp_ts / 1e4;
        dt          = median(diff(time_lfp_s));
        Fs_lfp      = 1/dt;
        trigOE.Fs_lfp      = Fs_lfp;
        trigOE.time_lfp_ts = time_lfp_ts;
        trigOE.time_lfp_s  = time_lfp_s;
        
        fprintf('  LFP Fs ~ %.2f Hz\n', Fs_lfp);
    else
        % sanity check same length
        if numel(t_ts) ~= numel(time_lfp_ts)
            warning('[extract_triggers_oe] time_lfp length differs for channel %d', ch);
        end
    end
    
    % TTL detection: reuse your detect_ttl_from_lfp_channel logic
    % v is a column vector
    [peak_idx, peak_val] = detect_ttl_from_lfp_channel(v);
    
    t_peak_ts = t_ts(peak_idx);
    t_peak_s  = t_peak_ts / 1e4;
    
    trigOE.(ttl_name).ch       = ch;
    trigOE.(ttl_name).peak_idx = peak_idx(:);
    trigOE.(ttl_name).peak_val = peak_val(:);
    trigOE.(ttl_name).t_raw_ts = t_peak_ts(:);
    trigOE.(ttl_name).t_raw_s  = t_peak_s(:);
    
    fprintf('  %s: detected %d TTL peaks on ch%d\n', ttl_name, numel(peak_idx), ch);
    
    % regularisation if we know expected count
    Nexp = [];
    if strcmp(ttl_name,'fus') && isfield(cfg,'N_fus_frames')
        Nexp = cfg.N_fus_frames;
    elseif strcmp(ttl_name,'face') && isfield(cfg,'N_face_frames')
        Nexp = cfg.N_face_frames;
    elseif strcmp(ttl_name,'eye') && isfield(cfg,'N_eye_frames')
        Nexp = cfg.N_eye_frames;
    end
    
    if ~isempty(Nexp) && ~isnan(Nexp) && Nexp > 0
        [t_reg_s, info] = regularize_ttl(t_peak_s, Nexp, ttl_name, cfg.plt);
        trigOE.(ttl_name).t_reg_s  = t_reg_s(:);
        trigOE.(ttl_name).reg_info = info;
        
        fprintf('    %s TTL regularized: raw=%d, final=%d, expected=%d\n', ...
            ttl_name, info.N_raw, info.N_final, info.N_expected);
    end
    
    if cfg.plt
        figure('Name', ['TTL ' ttl_name],'Color','w');
        subplot(2,1,1);
        plot(time_lfp_s, v); hold on;
        plot(t_peak_s, max(v)*ones(size(t_peak_s))*1.02,'.r');
        xlabel('Time (s)'); ylabel('LFP');
        title(sprintf('%s TTL channel %d', ttl_name, ch));
        subplot(2,1,2);
        histogram(diff(t_peak_s),50);
        xlabel('ITI (s)'); ylabel('count');
        title([ttl_name ' inter-trigger intervals']);
    end
end

end
