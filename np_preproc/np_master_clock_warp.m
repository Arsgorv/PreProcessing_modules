function [warpFn, info] = np_master_clock_warp(segFolder, npStreamRoot, masterStreamKey, opts)
% Build NP-clock -> master-clock time map for a segment.
%
% Strategy (tried in order):
%   0. If both streams' timestamps.npy already share the same clock
%      (Synchronizer was on at recording time), use timestamps.npy directly.
%   1. TTL event matching: pair rising-edge TTL events between NP and master
%      streams using IEI cross-correlation, build piecewise-linear interp.
%   2. Linear duration ratio fallback (current "Tier 1" behavior).
%
% Returns:
%   warpFn(t_np_s) -> t_master_s                  (vectorized function handle)
%   info.method, info.n_events, info.peak_corr, info.npDur_s, info.masterDur_s, ...

if nargin < 4, opts = struct(); end
if ~isfield(opts,'use_ttl_sync'),    opts.use_ttl_sync    = true; end
if ~isfield(opts,'sync_ttl_line'),   opts.sync_ttl_line   = []; end   % [] = use any rising edges
if ~isfield(opts,'min_ttl_events'),  opts.min_ttl_events  = 4; end
if ~isfield(opts,'min_ttl_corr'),    opts.min_ttl_corr    = 0.95; end

info = struct('method','none','n_events',0,'peak_corr',NaN, ...
              'npDur_s',NaN,'masterDur_s',NaN,'lag',NaN);

masterStreamRoot = oe_find_stream(segFolder, masterStreamKey);
if isempty(masterStreamRoot)
    warning('np_master_clock_warp: no master stream "%s" in %s; using identity', ...
        masterStreamKey, segFolder);
    warpFn = @(t) t;
    info.method = 'identity-no-master';
    return
end

[Fs_np,    nCh_np]    = oe_read_stream_info(npStreamRoot,     NaN, NaN);
[Fs_master, nCh_master] = oe_read_stream_info(masterStreamRoot, NaN, NaN);
if ~isfinite(Fs_np),    Fs_np    = 2500;  end
if ~isfinite(Fs_master), Fs_master = 30000; end

% raw stream sample counts
nNP     = oe_nSamples_stream(npStreamRoot,     nCh_np);
nMaster = oe_nSamples_stream(masterStreamRoot, nCh_master);
info.npDur_s     = double(nNP)     / Fs_np;
info.masterDur_s = double(nMaster) / Fs_master;

% ---------- 0) check timestamps.npy for pre-synced data ----------
ts_np_path     = fullfile(npStreamRoot,     'timestamps.npy');
ts_master_path = fullfile(masterStreamRoot, 'timestamps.npy');
if exist(ts_np_path,'file') == 2 && exist(ts_master_path,'file') == 2
    try
        ts_np_first = read_first_last_npy(ts_np_path);
        ts_m_first  = read_first_last_npy(ts_master_path);
        % If both first timestamps are >> 0 (master-clock anchored) AND
        % their last timestamps span ~the same window, the streams are
        % already synced. Heuristic: |npDur_via_ts - masterDur_via_ts| < 0.010 s.
        npDur_ts = ts_np_first(2) - ts_np_first(1);
        mDur_ts  = ts_m_first(2)  - ts_m_first(1);
        if abs(npDur_ts - mDur_ts) < 0.010
            % Use NP timestamps directly: NP-sample-index -> master-clock-time
            np_ts_all = read_npy(ts_np_path);   % seconds (master clock)
            np_idx_s  = (0:numel(np_ts_all)-1)' / Fs_np;  % NP-clock time
            warpFn = @(t) interp1(np_idx_s, double(np_ts_all), t, 'linear', 'extrap');
            info.method   = 'oe-synchronizer';
            info.n_events = numel(np_ts_all);
            return
        end
    catch
        % parsing failed; fall through to TTL/linear
    end
end

% ---------- 1) TTL event matching ----------
if opts.use_ttl_sync
    [np_evts_s,     np_lines]     = oe_read_ttl_events(npStreamRoot,     Fs_np);
    [master_evts_s, master_lines] = oe_read_ttl_events(masterStreamRoot, Fs_master);

    if ~isempty(np_evts_s) && ~isempty(master_evts_s)
        % Filter rising edges (state > 0). Optionally restrict to a single line.
        npRise = np_evts_s(np_lines > 0);
        mRise  = master_evts_s(master_lines > 0);
        if ~isempty(opts.sync_ttl_line)
            ln = opts.sync_ttl_line;
            npRise = np_evts_s(np_lines == ln);
            mRise  = master_evts_s(master_lines == ln);
        end

        if numel(npRise) >= opts.min_ttl_events && numel(mRise) >= opts.min_ttl_events
            [npPaired, mPaired, matchInfo] = match_event_trains(npRise, mRise);
            if numel(npPaired) >= opts.min_ttl_events && ...
               isfinite(matchInfo.peak_corr) && matchInfo.peak_corr >= opts.min_ttl_corr
                % piecewise linear with linear extrapolation
                warpFn = @(t) interp1(npPaired, mPaired, t, 'linear', 'extrap');
                info.method    = 'ttl';
                info.n_events  = numel(npPaired);
                info.peak_corr = matchInfo.peak_corr;
                info.lag       = matchInfo.lag;
                info.np_event_times_s     = npPaired(:);
                info.master_event_times_s = mPaired(:);
                return
            else
                warning('np_master_clock_warp: TTL match poor (n=%d, corr=%.3f); falling back', ...
                    numel(npPaired), matchInfo.peak_corr);
            end
        end
    end
end

% ---------- 2) linear duration-ratio fallback ----------
if info.npDur_s > 0
    ratio = info.masterDur_s / info.npDur_s;
else
    ratio = 1;
end
warpFn = @(t) t * ratio;
info.method = 'linear-ratio';
info.peak_corr = NaN;
info.n_events  = 0;
end

% =============================================================================
function [evt_s, lines] = oe_read_ttl_events(streamRoot, Fs)
% Read TTL rising/falling events for a stream. Returns event times (in
% stream-clock seconds) and line numbers (positive=rising, negative=falling).
evt_s = []; lines = [];

% events folder is sibling of "continuous": .../recordingX/events/<streamFolder>/TTL/
recRoot = fileparts(fileparts(streamRoot));         % .../recordingX/
[~, streamName, ext] = fileparts(streamRoot);
streamFolder = [streamName ext];
ttlPath = fullfile(recRoot,'events',streamFolder,'TTL');
if ~isfolder(ttlPath)
    % Try any TTL subfolder under events/<streamFolder>/
    base = fullfile(recRoot,'events',streamFolder);
    if isfolder(base)
        sub = dir(base);
        sub = sub([sub.isdir] & ~ismember({sub.name},{'.','..'}));
        for i = 1:numel(sub)
            cand = fullfile(sub(i).folder, sub(i).name);
            if exist(fullfile(cand,'sample_numbers.npy'),'file') == 2
                ttlPath = cand; break
            end
        end
    end
end
if ~isfolder(ttlPath), return; end

snFile = fullfile(ttlPath,'sample_numbers.npy');
stFile = fullfile(ttlPath,'states.npy');
if ~(exist(snFile,'file')==2 && exist(stFile,'file')==2), return; end

try
    sn = double(read_npy(snFile));
    st = double(read_npy(stFile));
catch
    return
end
if isempty(sn), return; end

% Sample numbers are stream-relative if events were written by OE 0.6+ on
% the same stream. Convert to seconds via stream Fs.
% Subtract first sample number of the continuous stream if present so that
% time 0 == continuous-data sample 0.
firstSn = read_first_continuous_sn(streamRoot);
if isfinite(firstSn)
    sn = sn - firstSn;
end
evt_s = sn / double(Fs);
lines = st;
end

function v = read_first_continuous_sn(streamRoot)
v = 0;
snFile = fullfile(streamRoot,'sample_numbers.npy');
if exist(snFile,'file') == 2
    try
        a = read_npy(snFile);
        v = double(a(1));
    catch
        v = 0;
    end
end
end

function v = read_first_last_npy(npyPath)
% Cheap helper: read full vector and pick first & last. NPY events are short.
a = read_npy(npyPath);
if isempty(a)
    v = [NaN NaN];
else
    v = double([a(1) a(end)]);
end
end

function [a_pair, b_pair, info] = match_event_trains(a, b)
% Match two rising-edge trains by inter-event-interval cross-correlation.
% a, b are sorted vectors of event times (any units, must be the same).
info = struct('lag',NaN,'peak_corr',NaN);
a_pair = []; b_pair = [];
if numel(a) < 4 || numel(b) < 4, return; end

a = a(:); b = b(:);
da = diff(a);
db = diff(b);
da = da - mean(da);
db = db - mean(db);

% Cross-correlate IEIs
m = min(numel(da), numel(db));
[c, lags] = xcorr(da, db, m-1, 'coeff');
[peak, kmax] = max(c);
lag = lags(kmax);

info.lag       = lag;
info.peak_corr = peak;

% Construct paired indices using the matched lag
if lag >= 0
    % a's IEI[i] matches b's IEI[i-lag]
    aIdx = (1+lag):min(numel(da), numel(db)+lag);
    bIdx = aIdx - lag;
else
    bIdx = (1-lag):min(numel(db), numel(da)-lag);
    aIdx = bIdx + lag;
end
if isempty(aIdx), return; end

% IEI index k corresponds to events [k, k+1] in the original train; take both.
a_evt_idx = unique([aIdx, aIdx(end)+1]);
b_evt_idx = unique([bIdx, bIdx(end)+1]);
a_evt_idx = a_evt_idx(a_evt_idx >= 1 & a_evt_idx <= numel(a));
b_evt_idx = b_evt_idx(b_evt_idx >= 1 & b_evt_idx <= numel(b));

n = min(numel(a_evt_idx), numel(b_evt_idx));
a_pair = a(a_evt_idx(1:n));
b_pair = b(b_evt_idx(1:n));
end