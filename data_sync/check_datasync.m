function check_datasync(datapath)

load(fullfile(datapath, 'Master_sync.mat'))
if contains(datapath, 'Edel')
    load(fullfile(datapath, 'ephys\LFPData\LFP26.mat'))
elseif contains(datapath, 'Kosichka')
    load(fullfile(datapath, 'ephys\SleepScoring_OBGamma.mat'), 'BrainPower')
end
fus_file = dir(fullfile(datapath, 'fUS\RP_data_*slice_A.mat'));
load(fullfile(datapath, 'fUS', fus_file.name))

figure
subplot(411)
% Epochs_it = {Epochs.trial_music ; Epochs.trial_ferret ; Epochs.stim_music; Epochs.stim_ferret};
% Epochs_it = {Epochs.trial_speech ; Epochs.trial_ferret ; Epochs.stim_speech; Epochs.stim_ferret};
Epochs_it = {Epochs.trial_speech ; Epochs.trial_music ; Epochs.stim_speech; Epochs.stim_music};
for i = 1:numel(Epochs_it)
    restr = Restrict(cat_tsd.data, Epochs_it{i});
    
    restr_data = Data(restr);
    restr_data = reshape(restr_data', cat_tsd.Nx, cat_tsd.Ny, size(restr_data, 1)');
    restr_data = squeeze(mean(restr_data, 1:2));
    
    hold on
    plot(Range(restr, 'min'), restr_data, '.');
    legend({'trial 1'; 'trial 2'; 'stim 1'; 'stim 2'})
end

subplot(412)
for i = 1:numel(Epochs_it)
    try
        restr = Restrict(LFP, Epochs_it{i});
    catch
        restr = Restrict( BrainPower.Power{1}, Epochs_it{i});
    end
    hold on
    plot(Range(restr, 'min'), Data(restr), '.');
    legend({'trial 1'; 'trial 2'; 'stim 1'; 'stim 2'})
end

subplot(413)
Epochs_fus = {Epochs.fus_preexp ; Epochs.fus_exp ; Epochs.fus_postexp};
for i = 1:numel(Epochs_fus)
    restr = Restrict(cat_tsd.data, Epochs_fus{i});
    
    restr_data = Data(restr);
    restr_data = reshape(restr_data', cat_tsd.Nx, cat_tsd.Ny, size(restr_data, 1)');
    restr_data = squeeze(mean(restr_data, 1:2));
    
    hold on
    plot(Range(restr, 'min'), restr_data, '.');
    legend({'pre exp'; 'exp'; 'post exp'})
end

subplot(414)
for i = 1:numel(Epochs_fus)
    try
        restr = Restrict(LFP, Epochs_fus{i});
    catch
        restr = Restrict( BrainPower.Power{1}, Epochs_fus{i});
    end
    hold on
    plot(Range(restr, 'min'), Data(restr), '.');
    plot(Range(BrainPower.Power{1}, 'min'), Data(BrainPower.Power{1})+1e4, '.');
    legend({'pre exp'; 'exp'; 'post exp'; 'raw'})
end

end
