function [binnedspk,chLST,clLST] = ARS_binspikes(spikeFolder,specDat,TORCmatfile,onsetFile)
% Parameters
load(TORCmatfile);
% load(onsetFile);
stim_start = onsetFile;
onsetTimes = reshape(stim_start,[30 106])';% Extract dimensions
num_repeats = size(data.stimuliNames, 1);  % Number of repetitions (106)
num_stimuli = size(data.stimuliNames, 2); % Number of stimuli (30)

% LOAD and BIN
fprintf(strcat(' load the files \n'))
SPKfnlst = dir([spikeFolder 'times*.mat']);
fileNb = length(SPKfnlst);
sr = 1/diff(specDat.time_vector(1:2));
bin_edges = [specDat.time_vector ; specDat.time_vector(end)+1/sr];
clear chLST clLST channelNumLst Y
n_cluster = 0;
clear binnedspk;
for FNum = 1:fileNb
    fileName = SPKfnlst(FNum).name;
    locusSTR = strfind(fileName,'times_C');
    channelNumLst(FNum) = str2num(fileName(locusSTR+7:end-4));
    load([spikeFolder fileName]);
    unique_clusters = unique(cluster_class(:,1)); % All unique cluster numbers
    unique_clusters(unique_clusters==0) = [];
    % Discard cluster 0 spikes
   
    for clNum = 1:length(unique_clusters)
        n_cluster = n_cluster + 1;
        cl_spikes = cluster_class(cluster_class(:,1)==clNum, 2)/1000;
        % Iterate over repeats
        for repeat = 1:num_repeats
            % Extract stimuli names and onset times for this repetition
            stim_names = data.stimuliNames(repeat, :);
            onset_times = onsetTimes(repeat, :);

            % Sort stimuli alphabetically within this repetition
            [sorted_names, sort_idx] = sort(stim_names);
            sorted_onset_times = onset_times(sort_idx);

            % Process each stimulus
            for stim_idx = 1:num_stimuli
                onset_time = sorted_onset_times(stim_idx);  % Sorted stimulus onset time
                stim_name = sorted_names{stim_idx};        % Sorted stimulus name
                window_start = onset_time;      % Start of binning window
                window_end = onset_time + bin_edges(end);       % End of binning window
                cl_spikesRS = cl_spikes(cl_spikes>window_start&cl_spikes<window_end);

                % Bin spikes within this window
                spike_counts = histcounts(cl_spikesRS-window_start, bin_edges);

                % Store binned spikes with format [time, stim, cluster]
                binnedspk(:,stim_idx,repeat,n_cluster) = spike_counts;
            end
        end

        if any(isnan(binnedspk(:,:,:,n_cluster)))
            disp(FNum)
        end
        chLST(n_cluster) = channelNumLst(FNum);
        clLST(n_cluster) = clNum;
    end
end
