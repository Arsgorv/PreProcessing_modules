clear
addpath(genpath('libs\zeta_test\'))

sub = 1;
subname = 'mozzarella';
conditions = {'train'};
% 
% 
% figure
% t = tiledlayout(3,1,"TileSpacing","compact");
for i = 1:length(conditions)
session = 'post';
condition = conditions{i};

%load(["D:\SpeechLearningProject\Electrophysiology\", session, '\', subname, '\sub', num2str(sub), '_300_8000Hz.mat'])
%load(["D:\SpeechLearningProject\Electrophysiology\CND", session, '\', subname, '\sub', num2str(sub), '_300_8000Hz.mat'])
load("D:\SpeechLearningProject\Electrophysiology\CND\dataSub_1.mat")
load("D:\SpeechLearningProject\Electrophysiology\post\mozzarella\sub1_300_8000Hz.mat")
load("D:\SpeechLearningProject\Electrophysiology\raw\mozzarella_CST_2024_11_29_1254_elmina.mat")
%load(['neural\',session,'\',subname,'_',condition,'\sub',num2str(sub),'_300_8000Hz.mat'])
%load(['dataCND\',session,'\',condition,'\dataSub_',num2str(sub),'_300_8000Hz.mat'])
%load(['dataCND\',session,'\',condition,'\dataStim_',condition,'.mat'])
stim = load("D:\SpeechLearningProject\StimGeneration\FinalAudio\stim_ref.mat")

presented_order = data.stimuliNames;
for i = 1:length(presented_order)
    temp_stim = presented_order{i};
    split_order = split(temp_stim, '_');
    id_order = split_order{2};
    split_order = split(id_order, '.');
    id_order = split_order{1};
    presented_order{i} = str2num(id_order);
end



% wavstimnames = cellfun(@(x) [x,'.wav'],stim.stimuli_names,'UniformOutput',false)';
% 
% [~,perm] = ismember(d.stimuli_names.stimuliNames,wavstimnames);
% 
% stim.stimuli_names = stim.stimuli_names(perm);
% stim.condIdxs = stim.condIdxs(perm);
% stim.stimIdxs = stim.stimIdxs(perm);
% stim.data = stim.data(:,perm);
% 
% onsets = stim.data(2,:);
%%
eventtimes = [];
trialEnds = 0;
for iTrial = 1:length(presented_order)
    trialStart = d.triggers(iTrial);
    if trialEnds > trialStart
        disp(iTrial)
        continue
    end
    onsettimes = (stim.(['stim_',num2str(presented_order{iTrial})]).duration(:, 1)*30)+1;
    eventtimes = [eventtimes; (trialStart+onsettimes)/30000];
    trialEnds = trialStart + max(onsettimes);
end

nChans = max(d.spikes(:,2));
%%
addpath(genpath("D:\SpeechLearningProject\Electrophysiology\libs\zetatest"))
sigs = zeros(1,nChans);
for i = 1:nChans
    channelIdxs = find(d.spikes(:,2)==i);
    spikeTimes = d.spikes(channelIdxs)/30000;
    [dblZetaP,sZETA,sRate,sLatencies] = zetatest(spikeTimes,eventtimes);
    disp(dblZetaP)
    sigs(i) = dblZetaP;
end
nexttile;
bar(sigs)
title(condition)

save(['D:\SpeechLearningProject\Electrophysiology\CND\sub',num2str(sub),'_',session,'_zetaP.mat'],"sigs")
end
xlabel(t,'Channel Number')
ylabel(t,'ZETA P val')
