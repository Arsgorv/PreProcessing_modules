function [peak_time_ts, info] = detect_ttl_from_lfp_channel(v, time_ts, lfp_chan)
% detect_ttl_from_lfp_channel
% Detect TTL pulses in an analog LFP trigger channel.
%
% INPUT
%   v         : [N x 1] analog samples (double or single)
%   time_ts   : [N x 1] time vector in ts units (1e-4 s), e.g. Range(LFP)
%   lfp_chan  : (optional) channel number for logging / warnings only
%
% OUTPUT
%   peak_time_ts : [nEvents x 1] timestamps (ts units, same as time_ts)
%   info         : struct with detection metadata

if nargin < 2
    error('detect_ttl_from_lfp_channel:NotEnoughInputs', ...
        'Need v and time_ts as inputs.');
end

if nargin < 3
    lfp_chan = NaN; % just for messages
end

v = v(:);
time_ts = time_ts(:);

if numel(v) ~= numel(time_ts)
    error('detect_ttl_from_lfp_channel:SizeMismatch', ...
        'v and time_ts must have the same length.');
end

if numel(v) < 3
    warning('detect_ttl_from_lfp_channel:ShortSignal', ...
        'LFP channel %g has too few samples to detect TTLs.', lfp_chan);
    peak_time_ts = [];
    info = struct('n_raw', numel(v), 'n_events', 0);
    return
end

dv = diff(v);
dv_max = max(dv);

if dv_max <= 0
    warning('detect_ttl_from_lfp_channel:FlatSignal', ...
        'LFP channel %g has no positive edges (max diff <= 0).', lfp_chan);
    peak_time_ts = [];
    info = struct('n_raw', numel(v), 'n_events', 0);
    return
end

% Threshold on derivative: large upward jump indicates TTL edge
thr_diff = dv_max / 3;
rise_idx = find(dv >= thr_diff);   % indices in dv (so sample index = rise_idx+1)

if isempty(rise_idx)
    warning('detect_ttl_from_lfp_channel:NoEdges', ...
        'No edges detected on LFP channel %g.', lfp_chan);
    peak_time_ts = [];
    info = struct('n_raw', numel(v), 'n_events', 0);
    return
end

% Group close indices into events
threshold_event = 4;        % minimum gap in samples to separate events
num_points = numel(rise_idx);

event_indices = 1;
for i = 2:num_points
    if (rise_idx(i) - rise_idx(i-1)) > threshold_event
        event_indices = [event_indices, i]; %#ok<AGROW>
    end
end
if event_indices(end) ~= num_points
    event_indices = [event_indices, num_points];
end

num_events   = numel(event_indices) - 1;
peak_indices = nan(num_events, 1);
peak_values  = nan(num_events, 1);

for i = 1:num_events
    idx_start = event_indices(i);
    if i < num_events
        idx_end = event_indices(i+1) - 1;
    else
        % last event goes until last detected rising edge
        idx_end = num_points;
    end
    
    seg_idx  = rise_idx(idx_start:idx_end);  % indices in dv
    % convert to indices in v (dv(k) ~ v(k+1)-v(k)), the jump is at k+1
    seg_idx_v = seg_idx + 1;
    seg_vals  = v(seg_idx_v);
    
    [mx, k] = max(seg_vals);
    peak_indices(i) = seg_idx_v(k);
    peak_values(i)  = mx;
end

% sanity
if any(diff(peak_indices) <= 0)
    warning('detect_ttl_from_lfp_channel:NonMonotonicPeaks', ...
        'Detected TTL peak indices are not strictly increasing.');
end

% Remove NaN events (just in case)
valid = ~isnan(peak_indices);
peak_indices = peak_indices(valid);
peak_values  = peak_values(valid);

% Optional: outlier removal on inter-event intervals
interpeak = diff(peak_indices);
out_idx   = [];

if ~isempty(interpeak)
    z      = zscore(double(interpeak));
    thr_z  = 3.5;
    out_idx = find(z > thr_z | z < -thr_z);
end

if ~isempty(out_idx)
    % simplest: drop the second event of each "weird" interval
    drop_idx = out_idx + 1;
    drop_idx(drop_idx > numel(peak_indices)) = [];
    peak_indices(drop_idx) = [];
    peak_values(drop_idx)  = [];
end

% map to time
peak_time_ts = time_ts(peak_indices);

info = struct();
info.n_raw        = numel(v);
info.n_events     = numel(peak_time_ts);
info.dv_max       = dv_max;
info.thr_diff     = thr_diff;
info.n_rise_idx   = numel(rise_idx);
info.n_outliers   = numel(out_idx);
info.peak_indices = peak_indices;
info.peak_values  = peak_values;

end
