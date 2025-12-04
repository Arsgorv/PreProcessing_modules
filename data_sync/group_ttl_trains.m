function [trial_start_ts,trial_stop_ts,info] = group_ttl_trains(pulse_ts,gap_sec)
% pulse_ts : column vector of all detected TTL pulses (ts units, 1e-4 s)
% gap_sec  : minimum gap between trials (seconds). Pulses separated by
%            more than this are considered different trials.

if nargin<2
    gap_sec = 0.2; % adjust if needed (0.5 s works well for Kosichka)
end

pulse_ts = sort(pulse_ts(:));

if isempty(pulse_ts)
    trial_start_ts = [];
    trial_stop_ts  = [];
    info = struct('n_pulses',0,'n_trials',0,'gap_ts',gap_sec*1e4);
    return
end

gap_ts = gap_sec*1e4;          % convert seconds ? ts units

% --- define trains: whenever the gap is large, a new trial starts
dt = diff(pulse_ts);
is_new_trial = [true; dt>gap_ts];  % IMPORTANT: keep first pulse
trial_id = cumsum(is_new_trial);

n_trials = trial_id(end);
trial_start_ts = zeros(n_trials,1);
trial_stop_ts  = zeros(n_trials,1);

for k = 1:n_trials
    idx = find(trial_id==k);
    trial_start_ts(k) = pulse_ts(idx(1));      % first pulse in train
    trial_stop_ts(k)  = pulse_ts(idx(end));    % last pulse in train
end

info = struct;
info.n_pulses  = numel(pulse_ts);
info.n_trials  = n_trials;
info.gap_ts    = gap_ts;
info.trial_id  = trial_id;
info.raw_ts    = pulse_ts;
info.stop_ts   = trial_stop_ts;
end
