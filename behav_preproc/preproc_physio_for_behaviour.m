function preproc_physio_for_behaviour(datapath, opts)
% preproc_physio_for_behaviour (React_Active)
%
% Creates: <datapath>\ephys\physio_for_behaviour.mat
% Resamples selected physio signals to DLC video time base (time_face).
%
% Saves (if available):
%   OB_gamma_rs, OB_gamma_log10_rs, OB_gamma_z_rs
%   RespPower_rs, RespPower_log10_rs, RespPhase_rs
%   EMGPower_rs, EMGPower_log10_rs, EMGPhase_rs
%   AccPower_rs, AccPower_log10_rs
%   HeartPower_rs, HeartPower_log10_rs, HeartPhase_rs   (optional, if heart LFP file exists)
%
% Requirements:
%   video/DLC_data.mat -> time_face
%   ephys/SleepScoring_OBGamma.mat -> BrainPower.Power{1} (tsd)

if nargin < 2, opts = struct(); end
if ~isfield(opts,'smootime_s'), opts.smootime_s = 0.006; end
if ~isfield(opts,'acc_file'), opts.acc_file = 'behavResources.mat'; end
if ~isfield(opts,'do_ob_fast_power'), opts.do_ob_fast_power = true; end
if ~isfield(opts,'ob_fast_smooth_s'), opts.ob_fast_smooth_s = 0.02; end % 20 ms
if ~isfield(opts,'ob_fast_gamma_band'), opts.ob_fast_gamma_band = [40 60]; end
if ~isfield(opts,'ob_fast_delta_band'), opts.ob_fast_delta_band = [0.5 4]; end

% Physiology
cfg = get_trigger_config(datapath);
if ~isfield(opts,'resp_lfp_file'), opts.resp_lfp_file = fullfile('LFPData',['LFP' num2str(cfg.respi) '.mat']); end
if ~isfield(opts,'emg_lfp_file'),  opts.emg_lfp_file  = fullfile('LFPData',['LFP' num2str(cfg.emg) '.mat']);  end
if ~isfield(opts,'ob_lfp_file'),   opts.ob_lfp_file   = fullfile('LFPData',['LFP' num2str(cfg.ob) '.mat']); end
if ~isfield(opts,'heart_lfp_file'), opts.heart_lfp_file = fullfile('LFPData',['LFP' num2str(cfg.heart) '.mat']); end 

% ------------------- video time base -------------------
V = load(fullfile(datapath,'video','DLC_data.mat'),'time_face');
% time_video_s = double(V.time_face(:)) / 1e4;
time_video_s = ra_time_to_seconds(double(V.time_face(:))); %fixed AG 15/02/2026

% Helper for runmean fallback
use_runmean = exist('runmean','file') == 2;

% Container to save
Sout = struct();

% ------------------- OB gamma/delta (BrainPower) -------------------
try
    S = load(fullfile(datapath,'ephys','SleepScoring_OBGamma.mat'),'BrainPower');
    BP = S.BrainPower;

    % find indices robustly
    iG = find(strcmpi(BP.signal_names, 'OB_gamma'), 1, 'first');
    iD = find(strcmpi(BP.signal_names, 'OB_delta'), 1, 'first');

    if isempty(iG)
        warning('OB_gamma not found in BrainPower.signal_names');
    else
        OB_gamma = BP.Power{iG};
        [t0,y0] = get_tsd_xy(OB_gamma);
        y_rs = interp1(t0, y0, time_video_s, 'linear', NaN);

        Sout.OB_gamma_rs = tsd(time_video_s*1e4, y_rs);
        Sout.OB_gamma_log10_rs = tsd(time_video_s*1e4, log10(max(y_rs, eps)));
        Sout.OB_gamma_z_rs = tsd(time_video_s*1e4, (y_rs - nanmean(y_rs)) / nanstd(y_rs));
    end

    if isempty(iD)
        warning('OB_delta not found in BrainPower.signal_names');
    else
        OB_delta = BP.Power{iD};
        [t0,y0] = get_tsd_xy(OB_delta);
        y_rs = interp1(t0, y0, time_video_s, 'linear', NaN);

        Sout.OB_delta_rs = tsd(time_video_s*1e4, y_rs);
        Sout.OB_delta_log10_rs = tsd(time_video_s*1e4, log10(max(y_rs, eps)));
        Sout.OB_delta_z_rs = tsd(time_video_s*1e4, (y_rs - nanmean(y_rs)) / nanstd(y_rs));
    end

catch ME
    warning('OB BrainPower preproc failed: %s', ME.message);
end

% ------------------- OB fast power (from OB LFP) -------------------
if isfield(opts,'do_ob_fast_power') && opts.do_ob_fast_power
    obPath = fullfile(datapath, 'ephys', opts.ob_lfp_file);
    if exist(obPath,'file')
        try
            A = load(obPath,'LFP'); OB_LFP = A.LFP;

            % gamma
            Fil = FilterLFP(OB_LFP, opts.ob_fast_gamma_band, 1024);
            [tF,xF] = get_tsd_xy(Fil);
            env = abs(hilbert(xF));

            dt = median(diff(tF));
            w = max(1, ceil(opts.ob_fast_smooth_s / dt));
            if use_runmean
                envS = runmean(env, w);
            else
                envS = movmean(env, w, 'omitnan');
            end

            G = tsd(tF*1e4, envS);
            G_rs = resample_tsd_to_video(G, time_video_s);
            Sout.OB_gamma_fast_rs = G_rs;
            Sout.OB_gamma_fast_log10_rs = tsd(time_video_s*1e4, log10(max(Data(G_rs), eps)));

            % delta
            Fil = FilterLFP(OB_LFP, opts.ob_fast_delta_band, 1024);
            [tF,xF] = get_tsd_xy(Fil);
            env = abs(hilbert(xF));

            dt = median(diff(tF));
            w = max(1, ceil(opts.ob_fast_smooth_s / dt));
            if use_runmean
                envS = runmean(env, w);
            else
                envS = movmean(env, w, 'omitnan');
            end

            Dp = tsd(tF*1e4, envS);
            Dp_rs = resample_tsd_to_video(Dp, time_video_s);
            Sout.OB_delta_fast_rs = Dp_rs;
            Sout.OB_delta_fast_log10_rs = tsd(time_video_s*1e4, log10(max(Data(Dp_rs), eps)));

        catch ME
            warning('OB fast power failed: %s', ME.message);
        end
    else
        warning('OB LFP file missing for fast power: %s', obPath);
    end
end

% ------------------- Respiration -------------------
respPath = fullfile(datapath, 'ephys', opts.resp_lfp_file);
if exist(respPath,'file')
    try
        A = load(respPath,'LFP'); respiration_LFP = A.LFP;

        Fil = FilterLFP(respiration_LFP,[0.1 1],1024);
        [tF,xF] = get_tsd_xy(Fil);

        dt = median(diff(tF));
        w = max(1, ceil(opts.smootime_s / dt));

        if use_runmean
            env = runmean(xF.^2, w);
        else
            env = movmean(xF.^2, w, 'omitnan');
        end

        RespPower = tsd(tF*1e4, env);
        RespPower_rs = resample_tsd_to_video(RespPower, time_video_s);

        Sout.RespPower_rs = RespPower_rs;
        Sout.RespPower_log10_rs = tsd(time_video_s*1e4, log10(max(Data(RespPower_rs), eps)));

        ph = angle(hilbert(zscore(xF))) * 180/pi + 180;
        RespPhase = tsd(tF*1e4, ph);
        Sout.RespPhase_rs = resample_tsd_to_video(RespPhase, time_video_s);
    catch ME
        warning('Resp preproc failed: %s', ME.message);
    end
else
    warning('Resp LFP file missing: %s', respPath);
end

% ------------------- EMG -------------------
emgPath = fullfile(datapath, 'ephys', opts.emg_lfp_file);
if exist(emgPath,'file')
    try
        A = load(emgPath,'LFP'); EMG_LFP = A.LFP;

        Fil = FilterLFP(EMG_LFP,[50 300],1024);
        [tF,xF] = get_tsd_xy(Fil);

        dt = median(diff(tF));
        w = max(1, ceil(opts.smootime_s / dt));

        if use_runmean
            env = runmean(xF.^2, w);
        else
            env = movmean(xF.^2, w, 'omitnan');
        end

        EMGPower = tsd(tF*1e4, env);
        EMGPower_rs = resample_tsd_to_video(EMGPower, time_video_s);

        Sout.EMGPower_rs = EMGPower_rs;
        Sout.EMGPower_log10_rs = tsd(time_video_s*1e4, log10(max(Data(EMGPower_rs), eps)));

        ph = angle(hilbert(zscore(xF))) * 180/pi + 180;
        EMGPhase = tsd(tF*1e4, ph);
        Sout.EMGPhase_rs = resample_tsd_to_video(EMGPhase, time_video_s);
    catch ME
        warning('EMG preproc failed: %s', ME.message);
    end
else
    warning('EMG LFP file missing: %s', emgPath);
end

% ------------------- Accelerometer (MovAcctsd) -------------------
accPath = fullfile(datapath, 'ephys',  opts.acc_file);
if exist(accPath,'file')
    try
        A = load(accPath,'MovAcctsd');
        MovAcctsd = A.MovAcctsd;

        [t0,y0] = get_tsd_xy(MovAcctsd);

        dt = median(diff(t0));
        w = max(1, ceil(opts.smootime_s / dt));

        if use_runmean
            env = runmean(y0.^2, w);
        else
            env = movmean(y0.^2, w, 'omitnan');
        end

        AccPower = tsd(t0*1e4, env);
        AccPower_rs = resample_tsd_to_video(AccPower, time_video_s);

        Sout.AccPower_rs = AccPower_rs;
        Sout.AccPower_log10_rs = tsd(time_video_s*1e4, log10(max(Data(AccPower_rs), eps)));
    catch ME
        warning('Accel preproc failed: %s', ME.message);
    end
else
    warning('behavResources missing (MovAcctsd): %s', accPath);
end

% ------------------- Heart (optional) -------------------
% You must set opts.heart_lfp_file to the correct LFPData\LFP?.mat that contains ECG/heart channel.
if ~isempty(opts.heart_lfp_file)
    heartPath = fullfile(datapath,'ephys', opts.heart_lfp_file);
    if exist(heartPath,'file')
        try
            A = load(heartPath,'LFP'); Heart_LFP = A.LFP;

            % default: broad band; edit as needed
            Fil = FilterLFP(Heart_LFP,[5 100],1024);
            [tF,xF] = get_tsd_xy(Fil);

            dt = median(diff(tF));
            w = max(1, ceil(opts.smootime_s / dt));

            if use_runmean
                env = runmean(xF.^2, w);
            else
                env = movmean(xF.^2, w, 'omitnan');
            end

            HeartPower = tsd(tF*1e4, env);
            HeartPower_rs = resample_tsd_to_video(HeartPower, time_video_s);

            Sout.HeartPower_rs = HeartPower_rs;
            Sout.HeartPower_log10_rs = tsd(time_video_s*1e4, log10(max(Data(HeartPower_rs), eps)));

            ph = angle(hilbert(zscore(xF))) * 180/pi + 180;
            HeartPhase = tsd(tF*1e4, ph);
            Sout.HeartPhase_rs = resample_tsd_to_video(HeartPhase, time_video_s);
        catch ME
            warning('Heart preproc failed: %s', ME.message);
        end
    else
        warning('Heart LFP file missing: %s', heartPath);
    end
end

% ------------------- save -------------------
outFile = fullfile(datapath,'ephys','physio_for_behaviour.mat');
save(outFile, '-struct', 'Sout', '-v7.3');

end

% ===================== minimal helpers ==================================
function [t_s, y] = get_tsd_xy(sig)
t_s = Range(sig,'s');
y = Data(sig);
t_s = double(t_s(:));
y = double(y(:));
end

function sig_rs = resample_tsd_to_video(sig, time_video_s)
t0 = Range(sig,'s');
y0 = Data(sig);
t0 = double(t0(:));
y0 = double(y0(:));
% y_rs = interp1(t0, y0, time_video_s, 'linear', 'extrap');
y_rs = interp1(t0, y0, time_video_s, 'linear', Nan); %fixed by AG 15/02/2026

sig_rs = tsd(time_video_s*1e4, y_rs);
end

function t_s = ra_time_to_seconds(t)
% ra_time_to_seconds
% Accepts either seconds or ts units (1e-4 s). Returns seconds.
t = t(:);
if numel(t) < 2
    t_s = t;
    return
end
dt = nanmedian(diff(t));
if ~isfinite(dt)
    t_s = t;
    return
end
% Heuristic: if dt > 1, it's almost certainly ts units (1e-4 s).
% (Typical DLC dt in seconds ~0.02; in ts ~200.)
if dt > 1
    t_s = t / 1e4;
else
    t_s = t;
end
end