function FilterAllIQ_IconeousOne
%{
    The purpose here is to filter my IQ data. The actual meat of the code
    was written by others; I am just wrapping it up in something that spits
    it out in a way useful for my filing system. 

    To understand the filtering, read Demene et al 2015. Basically we are
    doing an SVD on matrices of space vs time and removing the first 55 singular
    vectors. The rationale is that because tissue movements are more correlated
    spatially, those will be removed first. The method seems to work. 

   This is very slightly changed from my old version of the wrapper to accomidate
   the IconeousOne IQ data. 
%}
%% a temporary argument
  

globalArguments.ferretToAnalyze = 'Bichonnet';
globalArguments.dataHardDrive = 'NAS';

[globalArguments,dataAddresses] = loadWorkspaceModule(globalArguments);

filterDirectory = [dataAddresses.basefUSAddress globalArguments.ferretToAnalyze];

fUSdir = dir([filterDirectory]);
sliceFolders = fUSdir(~cellfun(@isempty,regexp({fUSdir.name},['\d\d$'])));
skipAlreadyDone = 0;
%sliceFolders = fUSdir(~cellfun(@isempty,regexp({fUSdir.name},'\d\d_'))); % for both ferrets
skipCounter = 0;

%% arguments related to choosing your filter
%   filterType = 'remove_large_vessels'; % types are: 'remove_tissue', 'keep_tissue', and 'remove_large_vessels'. 
  filterType = 'remove_tissue';
%   filterType = 'keep_all';
  switch filterType
    % filteredDataName must end with "Data.mat" for the code as written to work. Which should be fine. 
    case 'remove_tissue'
      filteredDataName = 'Filtered_Trial_Data.mat'; % as the default, it was named simply.
    case 'keep_tissue'
      filteredDataName = 'Tissue_Data.mat';
    case 'remove_large_vessels'
      filteredDataName = 'Capillary_Data.mat';
    case 'keep_all'
      filteredDataName = 'Unfiltered_Data.mat';
  end
%% arguments related to the parameters of the fourrier filter
  fourrierFilterType = 'butterworth'; %%% sixt order butterworth... 
  switch fourrierFilterType
    case 'celianStyle' 
      %{
        I intend for this not to be usable, but here for reference. 
        celian used a sixth-order butterworth with bandpass [20;60]...
        So from 20 to 60 Hz... He converted it to an IIR filter, which becomes
        of lower order, but which I disagree with because I will always prefer
        a filtfilt to avoid phase problems. I think he may have done the IIR
        to avoid similar problems by reducing the order but you should only do
        that if you apply the filter online... 
      
        Thinking more though... Because I am taking the power doppler, the arguments
        for use of filtfilt are not very strong because there is no need to maintain
        timing in the specified range. 
      
        That said I certainly need to verify the justification for these cutoffs... 
        
          -Specifically, 60 Hz. 
      %}
      Band = [20; 60];
      BP = num2str(Band);
      [z,p,k] = butter(6, Band'/(500/2));
      SOS = zp2sos(z,p,k);
      FilterString = ['Filter' BP(1,:) '-' BP(2,:)];
    case 'butterworth'
      Band = 60; %%% So as far as I can tell, the formula relating axial speed to cutoff is Cutoff (Hz) = (axial speed (mm/s))/20; This means that I don't want any low-pass applied at all if I can help it. 
      BP = num2str(Band);
      [z,p,k] = butter(6, Band'/(500/2),'low'); %%% 500 is fUS framerate. Was hardcoded by celian so probably won't change. But that's what it is. 
      SOS = zp2sos(z,p,k);
  end
%% use normcorre
  %{
    I must consider using normcore for each aquisition... Seems dangerous and
    maybe unnecessary but also it might be worthwhile. Not yet implemented though.
  %}
  use_normcorre = 0;
  warning('hack because of circumstance, to be deleted') % literally no clue what this refers to. Maybe I deleted it without deleting the warning?
for slice = 1:length(sliceFolders)
  %% count number of sessions, based on number of session folders
    sessionFolders = dir([sliceFolders(slice).folder '/' sliceFolders(slice).name]);
    sessionFolders = sessionFolders(~cellfun(@isempty,regexp({sessionFolders.name},'[sS]ession')));
    sessionFolders = sessionFolders(vertcat(sessionFolders.isdir)); % only want directories here
    numberOfSessions = length(sessionFolders); % number of fUS sessions available.
  for session = 1:numberOfSessions   
    %% count number of performanceBlocks in this session, based on number of scan directories...
      scanDirectories = dir([sessionFolders(session).folder '/' sessionFolders(session).name]);
      
      scanDirectories = scanDirectories(vertcat(scanDirectories.isdir));
      scanDirectories = scanDirectories(~cellfun(@isempty,regexp({scanDirectories.name},'_fus2D')));
      %scanDirectories = scanDirectories(~ismember({scanDirectories.name},{'.','..'})); % remove '.' and '..'
    for performanceBlock = 1:length(scanDirectories) % should scroll (in the 2022 case) active_fus2D, passive_fus2D, anatomy_fus2d
      workingDirfUS = dir([scanDirectories(performanceBlock).folder '/' scanDirectories(performanceBlock).name]);
      workingDirfUS = workingDirfUS(~cellfun(@isempty,regexp({workingDirfUS.name},'^mat_.*.mat$')));
      %% the IQ files should be in the format mat_#.at, which the number of digits varying. This unfortunately means they will be out of order. So we have to pull tricks.
        actualNumbers = zeros([length(workingDirfUS) 1]);
        for i = 1:length(workingDirfUS)
          
          actualNumbers(i) = str2double(workingDirfUS(i).name(5:(end-4)));
          ';';
        end
        [~,SortedIndexing] = sort(actualNumbers);
      workingDirfUS = workingDirfUS(SortedIndexing); % this works. 2 is 2, 15 is 15, 360 is 360...
      ';';
      if ~(skipAlreadyDone&exist([scanDirectories(performanceBlock).folder '/' scanDirectories(performanceBlock).name '_' filteredDataName],'file'))
        for i = 1:length(workingDirfUS)
          [workingDirfUS(i).folder '/' workingDirfUS(i).name]
                try
                  loadedData = load([workingDirfUS(i).folder '/' workingDirfUS(i).name]); 
                catch
                  warning('possibly a corrupted mat file, which might not be corrupted on your harddrive... If you replace it you can continue the code from here');
                  loadedData = load([workingDirfUS(i).folder '/' workingDirfUS(i).name]); 
                end
                if isfield(loadedData,'IQblock')
                  IQ = loadedData.IQblock;
                elseif isfield(loadedData,'IQ')
                  IQ = loadedData.IQ; %%% So the "IQ demodulation and beamforming steps" are done before saving the data I think...                   
                else
                  error('no applicable field name found within the data file, you should investigate');
                end
                    IQ_signal = IQ(:,:,1:min(500,size(IQ,3))); %%% seems to only take a maximum of 500 frames... Mine will have 300 or 200 always though. 
                    [nz, nx, nt] = size(IQ_signal);
                    IQ_signal = reshape(IQ_signal, [nz*nx, nt]);
                    cov_matrix = IQ_signal'*IQ_signal;
                    [Eig_vect, Eig_val]= eig(cov_matrix);
                    Eig_vect=fliplr(Eig_vect);
                    Eig_val=rot90(Eig_val,2);
                    M_ACP = IQ_signal*Eig_vect;    % on obtient les lambda*u
                    skipped_eig_val =[1:55];% 190:200];

                    IQF_tissu = M_ACP(:,skipped_eig_val)*Eig_vect(:,skipped_eig_val)';
                    IQF_tissu = reshape(IQF_tissu, [nz, nx, nt]);
                    IQ_signal = reshape(IQ_signal, [nz, nx, nt]);
                    switch filterType
                      case {'remove_tissue'}
                        IQF_corrected = IQ_signal-IQF_tissu;  
                      case 'remove_large_vessels'
                        IQF_corrected = IQ_signal-IQF_tissu;  
                        IQF_corrected = sosfilt(SOS, double(cat(3, IQF_corrected(:,:,100:-1:1), IQF_corrected)), 3);
                        IQF_corrected(:,:,1:100) = [];
                      case {'keep_tissue'}
                        IQF_corrected = IQF_tissu; %%% seems right to me...
                      case {'keep_all'}
                        IQF_corrected = IQ_signal; %%% no changes made. 
                    end
                    IQF = IQF_corrected;

                    %             if mod(ii, n_steps) == 1
                    %                 b0 = zeros(size(IQ_filtered,1), size(IQ_filtered,2), n_steps);
                    %             end
                    %             b0(:,:,1+mod((ii-1),n_steps)) = mean(abs(IQ_filtered).^2,3);
                    if i == 1
                        Doppler_image = zeros([size(mean(abs(IQF).^2,3)) length(workingDirfUS)]);
                    end
                    Doppler_image(:,:,i) = mean(abs(IQF).^2,3); %%% This is when you get the power doppler signal. "Known to be proportional to blood volume". I wonder if I can still use it for my tissue image, but we'll see. 
                    %figure;imagesc(Doppler_image(:,:,i))

           ';';

        end
        save([scanDirectories(performanceBlock).folder '/' scanDirectories(performanceBlock).name '_' filteredDataName],'Doppler_image')
        'lol';
      end
    end
  end
end
end

