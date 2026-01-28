function baphy_preprocessing(datapath)

%% CURRENTLY NOT USED VERSION OF THE SCRIPT: 20260123 AG

% BEWARE OF TWO ASSUMPTIONS:
% 1. Hits of NoMotor Target trials are forced to be NaN instead of 0 to avoid artificially decreasing hit rate
% 2. Hits of Early Target trials are forced to be 1 instead of 0 to avoid artificially decreasing hit rate

%% Initialization
fig_sanity = 0;
smoothWin = 5;    % set = 1 to turn smoothing off
minRunLen = 20;   % how many consecutive downward steps trigger the cut-off

if contains(datapath, 'Tvorozhok') || contains(datapath, 'Kosichka')
    trig_ch = 22;
else
    trig_ch = input('input trigger channel: ');
end

%% Load exptparams from the Mfile
mfilename = dir(fullfile(datapath, 'stim', '*.m'));
run(fullfile(datapath, 'stim', mfilename.name))

%% Parse events to find all TRIALSTART indices
nEvents = numel(exptevents);

trial_start_idx = find(contains({exptevents.Note}, 'TRIALSTART'));

Reference_idx = []; Target_idx = [];
stim_onset_times.Reference = []; stim_onset_times.Target    = [];
reward_onset_times.Target  = [];
nosound_idx  = []; nosound_type = {}; 

for ii = 1:numel(trial_start_idx)
    idx = trial_start_idx(ii);
    tr  = exptevents(idx).Trial;
    
    % Look ahead to find the next PreStimSilence event in this trial
    k = idx+1;  trial_type = ''; stimT = NaN; rewT = NaN;
%     trial_stim_time = NaN;
%     trial_reward_time = NaN;
%    
%     stim_found = false;
%     reward_found = false;
    while k <= nEvents && exptevents(k).Trial == tr
        note = exptevents(k).Note;
         % ---------- Reference / Target classification -------------------
        if contains(note, 'PreStimSilence')
            if contains(note, 'Reference'), trial_type = 'Reference'; Reference_idx(end+1) = tr;
            elseif contains(note, 'Target'),trial_type = 'Target';Target_idx(end+1) = tr;
            end
        end
        % ---------- stimulus onset --------------------------------------
        if isnan(stimT) && contains(note, 'Stim') && contains(note,'snippet_sequence')
            stimT = exptevents(k).StartTime;
        end
        % ---------- reward onset ----------------------------------------
        if isnan(rewT) && (contains(note, 'AUTOMATIC REWARD') ||...
                             contains(note, 'LICK,EARLY') ||...
                             contains(note, 'LICK,HIT,RIGHT'))
            rewT = exptevents(k).StartTime;
        end
        % ---------- nosound trial ---------------------------------------  
        if contains(note,'Stim') && isfield(exptevents(k),'Rove') && isequal(exptevents(k).Rove,0)
            nosound_idx(end+1) = tr;
            if contains(note,'Reference'), nosound_type{end+1}='Reference';
            else,                         nosound_type{end+1}='Target';
            end
        end
        k = k+1;
    end
    % store onsets -------------------------------------------------------
    switch trial_type
        case 'Reference'
            stim_onset_times.Reference(end+1) = stimT;
        case 'Target'
            stim_onset_times.Target(end+1)  = stimT;
            reward_onset_times.Target(end+1)= rewT;
    end
end

trial_id_struct.Reference = Reference_idx(:);
trial_id_struct.Target = Target_idx(:);
stim_onset_times.Target = stim_onset_times.Target(:);
stim_onset_times.Reference = stim_onset_times.Reference(:);
reward_onset_times.Target = reward_onset_times.Target(:);

%% just a sanity check to is if stim and reward are consistent across trials
if fig_sanity == 1
    stim_all = [stim_onset_times.Reference ; stim_onset_times.Target];
    reward_all = reward_onset_times.Target;
    stim_isi = diff(stim_all);
    reward_isi = diff(reward_all);
    
    figure;
    subplot(2,1,1);
    histogram(stim_isi,50);
    title('Stimulus Onset Interval');
    xlabel('Time (s)');
    ylabel('Count');
    subplot(2,1,2);
    histogram(reward_isi,50);
    title('Automatic Reward Interval');
    xlabel('Time (s)');
    ylabel('Count');
    disp(['Mean stim interval: ' num2str(mean(stim_isi)) ' s, Std: ' num2str(std(stim_isi))])
    disp(['Mean reward interval: ' num2str(mean(reward_isi)) ' s, Std: ' num2str(std(reward_isi))])
end

%% Load baphy triggers
load(fullfile(datapath, 'LFPData', ['LFP' num2str(trig_ch) '.mat']));
ttlIdx = thresholdIntervals(LFP, 2e4, 'Direction', 'Above');
ttlTime = Start(ttlIdx);

% quick check
% figure; subplot(312), plot(Range(LFP, 's'), Data(LFP)), hold on
% plot(ttlTime/1e4, 1e4, '*r');

% Compose the trial structure, mapping trial numbers to LFP trial onsets
trial_structure.trial_onset.Reference = ts(ttlTime(trial_id_struct.Reference));
trial_structure.trial_onset.Target = ts(ttlTime(trial_id_struct.Target));

% Stimulus onset
trial_structure.stim_onset.Reference = stim_onset_times.Reference;
trial_structure.stim_onset.Target = stim_onset_times.Target;

% Reward onset: only for target trials
trial_structure.reward_onset.Target = reward_onset_times.Target;

%% No-sound trial IDs
trial_id_struct.NoSound = unique(nosound_idx(:));
trial_id_struct.NoSoundType = nosound_type(1:3:end)';
trial_structure.trial_onset.NoSound = ts(ttlTime(unique(nosound_idx)));

%% No-motor trial IDs
if isfield(exptparams,'ProbeTrialLst')
    if isempty(exptparams.ProbeTrialLst)
        probeVec = [];
    else
        probeVec = str2double(strsplit(exptparams.ProbeTrialLst,','));
    end
    nomotor_idx = find(probeVec==1);
    exclusive_nomotor_idx = setdiff(nomotor_idx,nosound_idx);

    % label each exclusive nomotor trial as Ref / Tar
    nomotorType = cell(numel(exclusive_nomotor_idx),1);
    for kk = 1:numel(exclusive_nomotor_idx)
        tr = exclusive_nomotor_idx(kk);
        ix = find([exptevents.Trial]==tr & contains({exptevents.Note},'PreStimSilence'),1);
        if contains(exptevents(ix).Note,'Reference'), nomotorType{kk}='Reference';
        elseif contains(exptevents(ix).Note,'Target'), nomotorType{kk}='Target';
        else, nomotorType{kk}='Unknown';
        end
    end
    trial_id_struct.NomotorExclusive = exclusive_nomotor_idx(:);
    trial_id_struct.NomotorExclusiveType = nomotorType(:);
    trial_structure.trial_onset.NomotorExclusive = ts(ttlTime(exclusive_nomotor_idx));
end

%% Behaviour metrics (first lick, hit/lick rate, fatigue cut-off)

% Extract first licks
trial_structure.lick_onset.Target = exptparams.FirstLick.Tar(:);
trial_structure.lick_onset.Reference = exptparams.FirstLick.Ref(:);

% Extract lick rate
lickRate = [exptparams.Performance(1:end-1).LickRate].';
hitRate = [exptparams.Performance(1:end-1).HitRate].';

trial_structure.lick_rate.all = lickRate;
trial_structure.hit_rate.all = hitRate;
trial_structure.lick_rate.Reference = lickRate(trial_id_struct.Reference);
trial_structure.lick_rate.Target = lickRate(trial_id_struct.Target);
trial_structure.hit_rate.Reference = hitRate(trial_id_struct.Reference);
trial_structure.hit_rate.Target = hitRate(trial_id_struct.Target);

nTrials = numel(trial_structure.hit_rate.all);

% Extract hit data
temp_hit = nan(nTrials,1);
for i = 1:nTrials
    temp_hit(i) = exptparams.Performance(i).Hit;
end
temp_hit = temp_hit';
trial_structure.hit_d = temp_hit(trial_id_struct.Target);

%% Turn No-Motor Target to NaN
hit_raw = double(trial_structure.hit_d(:));   % 1 = hit, 0 = miss/other

if isfield(trial_id_struct,'NomotorExclusive') && ...
   isfield(trial_id_struct,'NomotorExclusiveType')

%     % indices of NomotorExclusive trials that are Target
    isNMTarget = strcmpi(trial_id_struct.NomotorExclusiveType,'Target');
%     idx_excl = trial_id_struct.NomotorExclusive(isNMTarget);
% 
%     % guard against out-of-range indices
%     idx_excl = idx_excl(idx_excl>=1 & idx_excl<=numel(hit_raw));
% 
%     hit_raw(idx_excl) = NaN;      % mark as “not applicable”
    nm_abs_idx = trial_id_struct.NomotorExclusive(isNMTarget); % absolute trial numbers
    % map absolute trial numbers to positions inside Target-only vector
    [inTar, nm_pos_in_tar] = ismember(nm_abs_idx, trial_id_struct.Target);
    nm_pos_in_tar = nm_pos_in_tar(inTar);
    if ~isempty(nm_pos_in_tar)
        trial_structure.hit_d(nm_pos_in_tar) = NaN;
    end
end

%% Turn Early Target hits to hits
outcomeList = {exptparams.Performance.ThisTrial}';
earlyGlobIdx = find(strcmpi(outcomeList,'Early'));      % indices in 1:nTrials

% Restrict to Target trials – hit_d has the same order as trial_id_struct.Target
[isEarlyInTar, locInTar] = ismember(earlyGlobIdx, trial_id_struct.Target);

if any(isEarlyInTar)
    tarPos = locInTar(isEarlyInTar);        % positions inside hit_d
    % sanity: make sure we don’t overwrite an existing hit
    alreadyHit = hit_raw(tarPos)==1;
    if any(alreadyHit)
        warning('Some ''Early'' trials already marked as Hit; leaving them as is.');
    end
    trial_structure.hit_d(tarPos) = 1;      % set early-lick trials to hit
end

%% Fatigue point
if smoothWin > 1
    hitRate = movmean(hitRate, smoothWin, 'omitnan');
end
% detect first minRunLen strictly decreasing run
dHr = diff(hitRate);        % negative if next trial is worse
isDown = dHr < 0;              % logical vector (length = nTrials-1)

% run-length encoding of isDown
switchRuns = [true; diff(isDown)~=0]; % starts of runs
runStart = find(switchRuns); % indices in isDown
runLen = diff([runStart; numel(isDown)+1]); % length of each run
runVal = isDown(runStart); % whether each run is “down”

% first run that is “down” AND long enough
ix = find(runVal & runLen >= minRunLen, 1, 'first');

%% define good trials
if isempty(ix)
    goodTrials = 1:numel(hitRate); % animal never phased out
    goodTrials = goodTrials';
else
    % runStart(ix) is an index in isDown (diff vector)
    firstBad = runStart(ix);         % this diff is the first negative of the run
    goodTrials = 1:firstBad;         % keep everything *before* the drop starts
    goodTrials = goodTrials';
end

trial_id_struct.goodTrials = goodTrials;

%% Find spout arrival tp
load(fullfile(datapath, 'video', 'DLC_data.mat'), 'spout_likelihood')
idxTarget = trial_id_struct.Target(:);
idxReference = trial_id_struct.Reference(:);
tAll(idxReference) = Range(trial_structure.trial_onset.Reference,'s');
tAll(idxTarget) = Range(trial_structure.trial_onset.Target,'s');
M = struct(); % container for all masks
M.mskTarget = false(nTrials,1);  M.mskTarget(idxTarget) = true;
M.mskRef = false(nTrials,1); M.mskRef(idxReference) = true;
idxRef = find(M.mskRef);
idxTar = find(M.mskTarget);

ref_onsets = tAll(idxRef);
tar_onsets = tAll(idxTar);

% ref_onsets = Range(trial_structure.trial_onset.Reference, 's');
% tar_onsets = Range(trial_structure.trial_onset.Target, 's');
n_ref = length(ref_onsets); 
n_tar = length(tar_onsets); 
% Define the properties of the spout arrival
thr = 0.6; 
minDur = 1*1e4; % if coord changes in less than this, then it's a glitch
mergeGap = 1*1e4; % fuse two if there is a short glitch

% Define the window of the spout arrival
arrLagNom = 3.95; % Lag with which the motor is supposed to start off from the trial onset
lagTol = 0.5; % Lag for spout movement, I wait spout in this time window
arrLagMin = arrLagNom - lagTol;
arrLagMax = arrLagNom + lagTol;

% Clean up the data a bit
temp = runmean(Data(spout_likelihood), 3); spout = tsd(Range(spout_likelihood), temp);
ivRaw = thresholdIntervals(spout, thr, 'Direction', 'Above');
ivClean = dropShortIntervals(ivRaw, minDur);
ivClean = mergeCloseIntervals(ivClean, mergeGap);

% Find spout_arrival times for each trial
stAll = Start(ivClean, 's');

spout_arrival_ref = nan(n_ref, 1);

for tr = 1:n_ref
   winBeg_ref = ref_onsets(tr) + arrLagMin;
   winEnd_ref = ref_onsets(tr) + arrLagMax;   
   
   cand = stAll(stAll >= winBeg_ref & stAll <= winEnd_ref);
   if ~isempty(cand)
       spout_arrival_ref(tr) = cand(1);
   end
end

spout_arrival_tar = nan(n_tar, 1);

for tr = 1:n_tar
   winBeg_tar = tar_onsets(tr) + arrLagMin;
   winEnd_tar = tar_onsets(tr) + arrLagMax;
   
   cand = stAll(stAll >= winBeg_tar & stAll <= winEnd_tar);
   if ~isempty(cand)
       spout_arrival_tar(tr) = cand(1);
   end
end

rel_spout_arrival = nan(nTrials, 1);
rel_spout_arrival(idxRef) = spout_arrival_ref - ref_onsets';
rel_spout_arrival(idxTar) = spout_arrival_tar - tar_onsets';
trial_structure.rel_spout_arrival(idxRef) = rel_spout_arrival(idxRef);
trial_structure.rel_spout_arrival(idxTar) = rel_spout_arrival(idxTar);
trial_structure.abs_spout_arrival_ref = spout_arrival_ref;
trial_structure.abs_spout_arrival_tar = spout_arrival_tar;

%% Find spout start
% check that spout arrival for ref is smaller than for tar
% find median ref and tar spout arrival, make a hist of how different these are across trials
% motor launch is at 3.95

%% Plot trials fig
idxTarget = trial_id_struct.Target(:);
idxReference = trial_id_struct.Reference(:);
idxGood = trial_id_struct.goodTrials(:);  
idxNoSound = trial_id_struct.NoSound(:);
nsType = trial_id_struct.NoSoundType;


% A. build onsets vector (seconds)
tAll = nan(nTrials,1);

tAll(idxReference) = Range(trial_structure.trial_onset.Reference,'s');
tAll(idxTarget) = Range(trial_structure.trial_onset.Target,'s');
tAll(idxNoSound) = Range(trial_structure.trial_onset.NoSound,'s');

% B. logical masks (1×nTrials)
mskTarget = false(nTrials,1);   mskTarget(idxTarget) = true;
mskReference = false(nTrials,1);   mskReference(idxReference)= true;
mskGood = false(nTrials,1);   mskGood(idxGood) = true;
mskNoSound = false(nTrials,1);   mskNoSound(idxNoSound) = true;

% C. sub-masks for No-Sound / No-Motor target-type 
% use the fact that nsType & nmType align 1-for-1 with their own idx lists
isNSTarget       = strcmp(nsType,'Target');
isNSRef          = strcmp(nsType,'Reference');

mskNoSound_Target = false(nTrials,1);  mskNoSound_Target(idxNoSound(isNSTarget)) = true;
mskNoSound_Reference = false(nTrials,1);  mskNoSound_Reference(idxNoSound(isNSRef)) = true;

if isfield(exptparams,'ProbeTrialLst')
    idxNoMotor = trial_id_struct.NomotorExclusive(:);
    nmType  = trial_id_struct.NomotorExclusiveType;
    isNMTarget = strcmp(nmType,'Target');
    isNMRef = strcmp(nmType,'Reference');
    
    tAll(idxNoMotor) = Range(trial_structure.trial_onset.NomotorExclusive,'s');
    mskNoMotor = false(nTrials,1);   mskNoMotor(idxNoMotor) = true;
    mskRegular = ~(mskNoSound | mskNoMotor);
    mskNoMotor_Target = false(nTrials,1);  mskNoMotor_Target(idxNoMotor(isNMTarget)) = true;
    mskNoMotor_Reference = false(nTrials,1);  mskNoMotor_Reference(idxNoMotor(isNMRef)) = true;
else
    mskRegular = ~(mskNoSound);
end

% C. visual palette
clr.Target = [0.98 0.50 0.45];       % salmon
clr.Reference = [0.00 0.39 0.00];       % dark-green

shp.Regular = '.'; shp.NoSound = 'x';   shp.NoMotor = '*';
sz.Good  = 10; sz.Other = 5;

% D. data vectors 
t  = tAll;                                 % (nTrials×1) seconds
hr = trial_structure.hit_rate.all(:);      % hit-rate

% E. plotting helper 
f1 = figure('Color', 'w','visible', 'off'); hold on; ylim([0 1]);
set(f1, 'Units', 'Pixels', 'Position', [1 1 1920 1080/3]);
plt = @(mask,shape,col,ms,label) ...
      plot(t(mask), hr(mask), shape,...
      'Color',col,...
      'MarkerSize',ms,...
      'DisplayName',label,...
      'LineStyle','none');

% G. 12 combinations
% regular
plt(mskRegular & mskTarget &  mskGood , shp.Regular, clr.Target, sz.Good , 'Target – good – regular');
plt(mskRegular & mskTarget & ~mskGood , shp.Regular, clr.Target, sz.Other, 'Target – other – regular');
plt(mskRegular & mskReference &  mskGood , shp.Regular, clr.Reference, sz.Good , 'Reference – good – regular');
plt(mskRegular & mskReference & ~mskGood , shp.Regular, clr.Reference, sz.Other, 'Reference – other – regular');

% no-sound
plt(mskNoSound_Target &  mskGood , shp.NoSound, clr.Target, sz.Good , 'Target – good – no-sound');
plt(mskNoSound_Target & ~mskGood , shp.NoSound, clr.Target, sz.Other, 'Target – other – no-sound');
plt(mskNoSound_Reference &  mskGood , shp.NoSound, clr.Reference, sz.Good , 'Reference – good – no-sound');
plt(mskNoSound_Reference & ~mskGood , shp.NoSound, clr.Reference, sz.Other, 'Reference – other – no-sound');

% no-motor
if isfield(exptparams, 'ProbeTrialLst')
    plt(mskNoMotor_Target &  mskGood , shp.NoMotor, clr.Target, sz.Good , 'Target – good – no-motor');
    plt(mskNoMotor_Target & ~mskGood , shp.NoMotor, clr.Target, sz.Other, 'Target – other – no-motor');
    plt(mskNoMotor_Reference & mskGood , shp.NoMotor, clr.Reference, sz.Good, 'Reference – good – no-motor');
    plt(mskNoMotor_Reference & ~mskGood , shp.NoMotor, clr.Reference, sz.Other, 'Reference – other – no-motor');
end
% H. cosmetics
ylabel('Hit rate'); xlabel('Time from session start (s)');
[~,session_name] = fileparts(datapath);
sgtitle(['Session ' strrep(session_name,'_','\_')]);

legend('Location','eastoutside', 'Interpreter','none');
grid on; box on;

%% Save to output file if desired
save(fullfile(datapath, 'stim', 'trial_structure.mat'), 'trial_structure', 'trial_id_struct');
saveas(f1, fullfile(datapath, 'stim', 'session_structure.png'))
disp('Trial parsing complete!');

end