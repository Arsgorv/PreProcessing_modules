function trigOE = extract_triggers_oe(datapath, cfg)
% extract_triggers_oe
% Detect TTL pulses for fUS, Baphy, cameras, etc from OpenEphys LFPData.
% Optionally regularize fUS and video TTLs using regularize_ttl.
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
%       plt          : 0/1  make quick sanity plots (default 0)
%
% OUTPUT
%   trigOE : struct with subfields:
%     .Fs_lfp           : sampling rate of LFP (Hz, estimated)
%     .time_lfp_ts      : time vector in ts units (1e-4 s)
%     .time_lfp_s       : time vector in seconds
%
%     For each TTL channel present (fus, baphy, face, eye, respi, heart, OneBox):
%        .<name>.ch        : channel index
%        .<name>.peak_idx  : sample indices of TTL peaks
%        .<name>.peak_val  : peak values
%        .<name>.t_raw_ts  : peak times in ts units (1e-4 s)
%        .<name>.t_raw_s   : peak times in seconds
%        .<name>.t_reg_s   : regularized times in seconds (if N_*_frames provided)
%        .<name>.reg_info  : struct returned by regularize_ttl

if nargin < 2 || isempty(cfg)
    error('extract_triggers_oe:NeedCfg', ...
        'Second argument cfg with channel indices is required.');
end

if ~isfield(cfg,'plt'), cfg.plt = 0; end

trigOE = struct;

% figure out animal name just for logging
animals = {'Edel','Chabichou','Kosichka','Tvorozhok','Shropshire','Brynza','Labneh','Mochi','Brayon'};
Session_params.animal = 'unknown';
for i = 1:numel(animals)
    if contains(datapath, animals{i})
        Session_params.animal = animals{i};
        break
    end
end

[~, sess_name] = fileparts(datapath);
disp(['[extract_triggers_oe] Session: ' sess_name ' (' Session_params.animal ')']);

ttl_names = {'fus','baphy','face_cam','eye_cam','respi','heart','OneBox'};
ch_fields = {'fus_ch','baphy_ch','face_ch','eye_ch','respi','heart','OneBox'};

Fs_lfp      = [];
time_lfp_ts = [];
time_lfp_s  = [];

for k = 1:numel(ttl_names)
    ttl_name = ttl_names{k};
    ch_field = ch_fields{k};

    if ~isfield(cfg, ch_field)
        continue
    end

    ch = cfg.(ch_field);
    if isempty(ch) || isnan(ch)
        continue
    end

    lfp_file = fullfile(datapath,'ephys','LFPData', ['LFP' num2str(ch) '.mat']);
    if ~exist(lfp_file, 'file')
        warning('[extract_triggers_oe] LFP file not found for %s channel %d: %s', ...
            ttl_name, ch, lfp_file);
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
        % define common time base from the first channel we load
        time_lfp_ts = t_ts(:);
        time_lfp_s  = time_lfp_ts / 1e4;
        dt          = median(diff(time_lfp_s));
        Fs_lfp      = 1/dt;

        trigOE.Fs_lfp      = Fs_lfp;
        trigOE.time_lfp_ts = time_lfp_ts;
        trigOE.time_lfp_s  = time_lfp_s;

        fprintf('  LFP Fs ~ %.2f Hz\n', Fs_lfp);
    else
        % sanity check: same length?
        if numel(t_ts) ~= numel(time_lfp_ts)
            warning('[extract_triggers_oe] time_lfp length differs for channel %d', ch);
        end
    end

    % TTL detection: use detect_ttl_from_lfp_channel
    [peak_time_ts, info_det] = detect_ttl_from_lfp_channel(v, t_ts, ch, Inf);

    if contains(ttl_name,'baphy') && contains(datapath, 'Kosichka')
        disp('Grouping TTL trains')
        gap_sec = 0.1; 
        % compress trains ? one start and one stop per trial
        [trial_start_ts,trial_stop_ts,info_group] = group_ttl_trains(peak_time_ts);
        
        % for Baphy alignment, use only trial starts
        peak_time_ts = trial_start_ts;
        
        trigOE.(ttl_name).n_raw_pulses   = info_group.n_pulses;
        trigOE.(ttl_name).n_trials       = info_group.n_trials;
        trigOE.(ttl_name).trial_start_ts = trial_start_ts;
        trigOE.(ttl_name).trial_stop_ts  = trial_stop_ts;
        trigOE.(ttl_name).gap_ts         = info_group.gap_ts;
    end
    peak_idx = info_det.peak_indices(:);
    peak_val = info_det.peak_values(:);

    t_peak_ts = peak_time_ts(:);
    t_peak_s  = t_peak_ts / 1e4;

    trigOE.(ttl_name).ch       = ch;
    trigOE.(ttl_name).peak_idx = peak_idx;
    trigOE.(ttl_name).peak_val = peak_val;
    trigOE.(ttl_name).t_raw_ts = t_peak_ts;
    trigOE.(ttl_name).t_raw_s  = t_peak_s;
%     trigOE.(ttl_name).n_frames  = info_det.n_events;
%     trigOE.(ttl_name).n_trigs  = numel(trigOE.(ttl_name).t_raw_s(:));

    fprintf('  %s: detected %d TTL peaks on ch%d\n', ttl_name, numel(peak_idx), ch);

    if cfg.plt
        figure('Name', ['TTL ' ttl_name],'Color','w');
        subplot(2,1,1);
        plot(time_lfp_s, v);
        hold on;
        plot(t_peak_s, max(v)*ones(size(t_peak_s))*1.02,'.r');
        xlabel('Time (s)');
        ylabel('LFP');
        title(sprintf('%s TTL channel %d', ttl_name, ch));

        subplot(2,1,2);
        histogram(diff(t_peak_s),50);
        xlabel('ITI (s)');
        ylabel('count');
        title([ttl_name ' inter-trigger intervals']);
    end
end

end