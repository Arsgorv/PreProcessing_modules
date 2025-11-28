function run_normcorre_batch
%{
    "run_normcorre_batch" is my function for implementing the NoRMCorre template
    matching on fUS data...

    A notable aspect of this code is that I have allowed for the splitting of the
    dataset for which you calculate the nonrigid shifts from the dataset to which
    you apply those shifts. My main application for this has been "calculate shifts
    while I can see the large blood vessels, and then apply those shifts to the
    data where I filtered those out". I intend to also someday calculate the shifts
    with an image of the tissue, then apply the shifts to the vessel image...


    "run_normcorre_batch" makes use of many functions in the NoRMCorre package.
    Here is the citation:
      Eftychios A. Pnevmatikakis and Andrea Giovannucci, *NoRMCorre: An online algorithm for piecewise rigid motion correction of calcium imaging data*, Journal of Neuroscience Methods, vol. 291, pp 83-94, 2017; doi: [https://doi.org/10.1016/j.jneumeth.2017.07.031](https://doi.org/10.1016/j.jneumeth.2017.07.031)

    
%}


%% debug variables
dontsaveanything = 0; % turn this on if you are debugging and don't want to save.
doMovie = 1; % decides whether to make a movie of the template-match... Takes time, so turn off if you don't need it.
loadDataForMovie = 0; % if normcorre data is available, load it instead of recalculating template-match for the movie.
showMovie = 0; % decides whether to show the movie. Often it is better to wait and load the avi, but if you want to see something quickly...
dontSaveResultsAfterMovie = 0; % allows you to save the movie without saving normcorre (I think). It's fine to keep at 0 anyway.
if showMovie
    visibleFigures = 'on';
else
    visibleFigures = 'off';
end
dontactuallyrunnormcore = 0; % if 1, I reshape the data into the format as if normcorre had been run, but don't run it. Useful for datasets where bulging isn't a problem...
if dontactuallyrunnormcore
    warning('not actually running normcorre. Might be intentional')
end
saveCroppedDimensions = 0; %%% set to 1 if you want to make the dimension croppings for each image.
makeNewCroppings = 0; %%% set to 1 to make new croppings even if you already have some saved.

%% load basic addresses into workspace and add the address of the CCAcorrectionfUS function
%   ferretToAnalyze = 'Laguiole';
%   dataHardDrive = 'WORKDRIVE2';
%    ferretToAnalyze = 'Laguiole';
%    dataHardDrive = 'Elements';
%globalArguments.ferretToAnalyze = 'PigouilleOBG';
%globalArguments.ferretToAnalyze = 'Munster';
%globalArguments.dataHardDrive = 'OBGDrive';
globalArguments.ferretToAnalyze = 'Edel';
globalArguments.dataHardDrive = 'NAS';

[globalArguments,dataAddresses] = loadWorkspaceModule(globalArguments);
addpath(genpath([dataAddresses.homePath 'PHD_Code/DownloadedCode/NoRMCorre-master']))

%% normcorre variables which are not down in loadWorkspaceModule
globalArguments = normcorre_set_arguments(globalArguments);
for slice = 1:length(dataAddresses.sliceFolders) %%% selecting a specific slice to show off.
    %% count number of sessions, based on number of session folders
    sessionFolders = dir([dataAddresses.BehaveAddress globalArguments.ferretToAnalyze '/' dataAddresses.sliceFolders(slice).name]);
    sessionFolders = sessionFolders(~cellfun(@isempty,regexp({sessionFolders.name},'[sS]ession')));
    sessionFolders = sessionFolders(vertcat(sessionFolders.isdir)); % only want directories here
    numberOfSessions = length(sessionFolders); % number of fUS sessions available.
    for session = 1:numberOfSessions
        
        %% count number of performanceBlocks in this session, based on number of trialInfo files in
        trialInfoFiles = dir([sessionFolders(session).folder '/' sessionFolders(session).name]);
        trialInfoFiles = trialInfoFiles(~cellfun(@isempty,regexp({trialInfoFiles.name},'trialinfo')));
        trialInfoFileNumber = length(trialInfoFiles); % number of fUS sessions available.
        %% load file containing performanceBlock labels, which are unique among sessions for each slice
        load([trialInfoFiles(1).folder '/specificPerformanceBlockStrings.mat'],'specificPerformanceBlockStrings');
        for performanceBlock = 1:trialInfoFileNumber
            
            %% setup variablesToKeep
            variablesToKeep = [];
            variablesToKeep = who;
            %% get the trialInfoFile loaded and confirm current performanceBlock type.
            load([trialInfoFiles(performanceBlock).folder '/' trialInfoFiles(performanceBlock).name]); % load trialinfo file
            
            currentPerformanceBlockString = specificPerformanceBlockStrings{performanceBlock};
            sessionNameString = sessionFolders(session).name;
            
            %% find all fUS trial data and make sure you only keep those you have PerformanceData for.
            dataAddresses = define_workingDirfUS(currentPerformanceBlockString,sessionNameString,slice,globalArguments,dataAddresses); %%% defines "workingDirfUS" which will contain the folders with each trials' fUS data.
            switch globalArguments.ferretToAnalyze
                case {'Munster','Cabecou','Banon','Laguiole'} % these are the "old" ferrets.
                    trialFolders = dataAddresses.workingDirfUS(~cellfun(@isempty,regexp({workingDirfUS.name},'\d\d\d$')));
                    trialsWithfUSData = 2:length(trialFolders); % Here I deal with the first trial lack-of-data problem
                    if max(trialsWithfUSData) < size(PerformanceData,1)
                        reducedPerformanceData = PerformanceData(trialsWithfUSData,:);
                    elseif max(trialsWithfUSData) >= size(PerformanceData,1)
                        reducedPerformanceData = PerformanceData(2:end,:);
                    else
                        warning('if this happens you probably have an empty trialsWithfUSData and should rename some folders')
                        warning(['Specifically: check slice ' sliceFolders(slice).name ' ' sessionNameString ' ' currentPerformanceBlockString])
                        continue
                    end
                    %% gather all trials together
                    usableTrials = [0;zeros([(min(length(trialsWithfUSData),size(reducedPerformanceData,1))) 1])]; % first trial always 0;
                    for trialCount = 1:(min(length(trialsWithfUSData),size(reducedPerformanceData,1))) % we want to iterate through every trial and load up the fUS data. It will be worthwhile to keep track of n_repetitions and n_stim...
                        if trialCount == 1
                            b0 = [];
                            b02 = [];
                        end
                        trial = trialsWithfUSData(trialCount);
                        %% load b0, the image to build template off of
                        tempTrialData = load([trialFolders(trial).folder ['/'] trialFolders(trial).name '/' normCorreDataName]); % load the filtered Trial Data. It includes the aquisition times and the fUS data, but I will only be grabbing the fUS data
                        tempTrialData = tempTrialData.Doppler_image;
                        if size(tempTrialData,3) >= minFramesForGoodRecording
                            tempTrialData = tempTrialData(:,:,framesIncluded);
                            %% initialize arrays if this is the first loop
                            
                            b0 = cat(4,b0,tempTrialData);
                            usableTrials(trial) = 1;
                        else
                            usableTrials(trial) = 0;
                        end
                        %% load b02, the image to apply shifts to.
                        tempTrialData = load([trialFolders(trial).folder ['/'] trialFolders(trial).name '/' applyDataName]); % load the filtered Trial Data. It includes the aquisition times and the fUS data, but I will only be grabbing the fUS data
                        tempTrialData = tempTrialData.Doppler_image;
                        if size(tempTrialData,3) >= minFramesForGoodRecording
                            tempTrialData = tempTrialData(:,:,framesIncluded);
                            %% initialize arrays if this is the first loop
                            
                            b02 = cat(4,b02,tempTrialData);
                        end
                        
                        
                        
                    end
                case {'Bichonnet','Pigouille','PigouilleOBG'}
                    
                    if exist([dataAddresses.basefUSAddressWORKDRIVE globalArguments.ferretToAnalyze '/' dataAddresses.sliceFolders(slice).name '/' sessionNameString '/' ],'dir')
                        %% first, load up the template data
                        switch globalArguments.normcorreFilterType
                            case 'defaultSVD'
                                if exist([dataAddresses.workingDirfUS(1).folder '/' [specificPerformanceBlockStrings{performanceBlock} '_fus2D.source.scan']],'file')
                                    nfo = h5info([dataAddresses.workingDirfUS(1).folder '/' [specificPerformanceBlockStrings{performanceBlock} '_fus2D.source.scan']]); % Get HDF5 info, including data set names
                                    
                                    b0 = squeeze(h5read(nfo.Filename,'/Data')); % Load raw data
                                    b0 = permute(b0,[2 1 3]);
                                    anatomyImage_allFrames = mean(b0,3);
                                else
                                    warning(['missing file ' [dataAddresses.workingDirfUS(1).folder '/' [specificPerformanceBlockStrings{performanceBlock} '_fus2D.source.scan']] ', but maybe it isn''t supposed to exist']);
                                    continue
                                end
                            otherwise % all other filter types can be handled as below:
                                tempLoadedFile = load(([dataAddresses.workingDirfUS(1).folder '/' [specificPerformanceBlockStrings{performanceBlock} '_fus2D_' dataAddresses.applyDataSaveFolderName '.mat']]));
                                b0 = tempLoadedFile.Doppler_image;
                                clear tempLoadedFile;
                        end
                        %% then, load up the data which will be used.
                        switch globalArguments.dataFilterType
                            case globalArguments.normcorreFilterType % if same
                                b02 = b0;
                            case 'defaultSVD'
                                if exist([dataAddresses.workingDirfUS(1).folder '/' [specificPerformanceBlockStrings{performanceBlock} '_fus2D.source.scan']],'file')
                                    nfo = h5info([dataAddresses.workingDirfUS(1).folder '/' [specificPerformanceBlockStrings{performanceBlock} '_fus2D.source.scan']]); % Get HDF5 info, including data set names
                                    
                                    b0 = squeeze(h5read(nfo.Filename,'/Data')); % Load raw data
                                    b0 = permute(b0,[2 1 3]);
                                    anatomyImage_allFrames = mean(b0,3);
                                else
                                    warning(['missing file ' [dataAddresses.workingDirfUS(1).folder '/' [specificPerformanceBlockStrings{performanceBlock} '_fus2D.source.scan']] ', but maybe it isn''t supposed to exist']);
                                    continue
                                end
                            otherwise % all other filter types can be handled as below:
                                if exist([dataAddresses.workingDirfUS(1).folder '/' [specificPerformanceBlockStrings{performanceBlock} '_fus2D_' dataAddresses.applyDataSaveFolderName '.mat']],'file')
                                    tempLoadedFile = load(([dataAddresses.workingDirfUS(1).folder '/' [specificPerformanceBlockStrings{performanceBlock} '_fus2D_' dataAddresses.applyDataSaveFolderName '.mat']]));
                                    b02 = tempLoadedFile.Doppler_image;
                                    clear tempLoadedFile;
                                    'bluhahih';
                                else
                                    disp(['no file corresponding to this filter for ' [dataAddresses.workingDirfUS(1).folder '/' [specificPerformanceBlockStrings{performanceBlock} '_fus2D_' dataAddresses.applyDataSaveFolderName '.mat']]]);
                                    continue
                                end
                        end
                        
                    else
                        continue
                    end
                    ';';
            end
            
            %% ascertain that you are inputting a cropped image. Don't want to corrupt the out-region, and normcorre isn't designed for it. Specifically, because I have such large overlaps horizontally, I will corrupt the lateral edges if I include the under-cement here.
            if exist([trialInfoFiles(1).folder '/pre_normcorre_cropDimensions.mat' ],'file')
                load([trialInfoFiles(1).folder '/pre_normcorre_cropDimensions.mat' ],'pre_normcorre_cropDimensions')
            else
                warning('missing cropped file, making new one. Cancel out of code if you don''t want this')
                makeNewCroppings = 1;
            end
            if makeNewCroppings
                figure;imagesc(b0(:,:,1,1));
                good_crop = 0;
                while ~good_crop
                    left_adjust = input('Enter left adjust: ');
                    right_adjust = input('Enter right adjust: ');
                    bottom_adjust = input('Enter bottom adjust: ');
                    if isempty(left_adjust)
                        left_adjust = 10;
                    end
                    if isempty(right_adjust)
                        right_adjust = 30;
                    end
                    if isempty(bottom_adjust)
                        bottom_adjust = 20;
                    end
                    figure;imagesc(b0(1:(size(b0,1) - bottom_adjust),left_adjust:(size(b0,2)-right_adjust),1,1));
                    good_crop = input('Good crop? Enter 1 for yes, 0 or nothing for no: ');
                    if isempty(good_crop)
                        good_crop = 0;
                    end
                end
                pre_normcorre_cropDimensions = {[1:(size(b0,1) - bottom_adjust)]',[left_adjust:(size(b0,2)-right_adjust)]'};
                save([trialInfoFiles(1).folder '/pre_normcorre_cropDimensions.mat' ],'pre_normcorre_cropDimensions')
                if ~isempty(variablesToKeep)
                    clearvars('-except',variablesToKeep{:})
                end
                
                continue
                
            end
            if saveCroppedDimensions % just saving cropped dimensions.
                if ~isempty(variablesToKeep)
                    clearvars('-except',variablesToKeep{:})
                end
                continue
            end
            %% set up save address
            saveAddress = [dataAddresses.basefUSAddress globalArguments.ferretToAnalyze '/' globalArguments.ferretToAnalyze '_RepositionedData/' dataAddresses.normcorre_folders '/'];
            
            
            switch globalArguments.ferretToAnalyze
                case {'Munster','Cabecou','Banon','Laguiole'}
                    saveString = ['Repositioned_' dataAddresses.sliceFolders(slice).name '_' currentPerformanceBlockString];
                case {'Bichonnet','Pigouille','PigouilleOBG'}
                    % need session names for these guys...
                    saveString = ['Repositioned_' dataAddresses.sliceFolders(slice).name '_' sessionNameString '_' currentPerformanceBlockString];
                    
            end
            saveName = [saveString '.mat'];
            demoFigVideoAddress = [saveAddress '/' saveString '_demoFig.avi'];
            if loadDataForMovie&doMovie&exist([saveAddress saveName],'file')
                load([saveAddress saveName],'Y_in','Y_out','M2','shifts2','template2','options_nonrigid','pre_normcorre_cropDimensions','usableTrials','b0_in_initialSize','framesIncluded');
            else
                %% define b0_in
                b0_in = b0(pre_normcorre_cropDimensions{1},pre_normcorre_cropDimensions{2},:,:);
                b02_in = b02(pre_normcorre_cropDimensions{1},pre_normcorre_cropDimensions{2},:,:);
                %% set loop-variable arguments
                if globalArguments.normcorreArguments.vary_grid_size
                    switch globalArguments.normcorreFilterType
                        case {'tissue'}
                            globalArguments.normcorreArguments.grid_size = [round(size(b0_in,1)/10),6];
                        case 'defaultSVD'
                            globalArguments.normcorreArguments.grid_size = [size(b0_in,1),6];
                    end
                end
                %% rehape b0 for use with normcorre
                b0_in_initialSize = size(b0_in);
                Y_in = reshape(b0_in,[size(b0_in,[1 2]) prod(size(b0_in,[3 4]))]);
                Y_out = reshape(b02_in,[size(b02_in,[1 2]) prod(size(b02_in,[3 4]))]);
                switch globalArguments.normcorreFilterType
                    %% something I'm trying here: transforming the data before template matching to modulate prioritization. Want to reduce dynamic range...
                    case {'tissue','keep_all'}
                        Y_in = sqrt(sqrt(sqrt(Y_in)));
                end
                switch globalArguments.dataFilterType
                    %% when I want to test success of the above switch statement, I want to do the same to Y_out
                    case {'tissue','keep_all'}
                        Y_out = sqrt(sqrt(sqrt(Y_out)));
                end
                %% apply arguments via NoRMCorreSetParms
                options_nonrigid = NoRMCorreSetParms('d1',size(Y_in,1),'d2',size(Y_in,2),'grid_size',globalArguments.normcorreArguments.grid_size,'mot_uf',globalArguments.normcorreArguments.mot_uf,'max_shift',globalArguments.normcorreArguments.max_shift,'max_dev',globalArguments.normcorreArguments.max_dev,'init_batch',globalArguments.normcorreArguments.init_batch,'overlap_pre',globalArguments.normcorreArguments.overlap_pre,'overlap_post',globalArguments.normcorreArguments.overlap_post,'shifts_method',globalArguments.normcorreArguments.shifts_method,'min_diff',globalArguments.normcorreArguments.min_diff,'correct_bidir',globalArguments.normcorreArguments.correct_bidir);
                if ~dontactuallyrunnormcore
                    %% run run_normcorre_batch
                    [M2,shifts2,template2,options_nonrigid] = normcorre_batch_twoArrays(Y_in,Y_out,options_nonrigid);
                else
                    M2 = Y_out;
                    shifts2 = 0; % not the right format.. But this way I can save it, and also, it correctly represented that nothing is shifted.
                    template2 = nanmean(M2(:,:,1:2),3);
                end
            end
            %% stuff from the normcorre demo, to evaluate successes.
            %% plot a movie with the results
            if doMovie
                %% compute metrics
                %                   nnY = quantile(Y_out(:),0.005);
                %                   mmY = quantile(Y_out(:),0.995);
                nnY = quantile(Y_out(:),0.1);
                mmY = quantile(Y_out(:),0.9);
                
                %[cY,mY,vY] = motion_metrics(Y_in,10);
                [cM2,mM2,vM2] = motion_metrics(M2,10);
                T = length(cM2);
                
                demoFig = figure('visible',visibleFigures);
                set(demoFig,'Position',[50 50 1200 800]);
                makeF = 1;
                for t = 1:1:T
                    %shifts2(t).shifts;
                    
                    subplot(131);imagesc(sqrt(Y_out(:,:,t)),sqrt([nnY,mmY])); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
                    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('hot');
                    subplot(132);quiv = quiver(flipud(shifts2(t).shifts_up(:,:,1,2)),flipud(shifts2(t).shifts_up(:,:,1,1)),'AutoScale','off');
                    %ylim(quiv.Parent,[-(max_shift(1)+1) (max_shift(1)+1)]); %%% do this if you want to focus on the y axis.
                    ylim(quiv.Parent,[-10 10]);
                    axis equal; %%% do this if you want the directions to be correct. I only care about the y axis, so no need.
                    xlabel('shifts_up','fontsize',14,'fontweight','bold');
                    
                    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('hot');
                    subplot(133);imagesc(sqrt(M2(:,:,t)),sqrt([nnY,mmY])); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
                    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('hot');
                    set(gca,'XTick',[],'YTick',[]);
                    drawnow;
                    temp =  getframe(demoFig);
                    F(t) = temp;
                    if makeF
                        makeF = 0;
                        F(T+1) = temp;
                        F(T+1) = [];
                    end
                end
                
                writerObj = VideoWriter(demoFigVideoAddress);
                set(writerObj,'FrameRate',30); % this is far faster than the actual session... See if it is tolerable.
                open(writerObj);
                writeVideo(writerObj, F)
                close(writerObj);
            end
            
            %% save the results
            if (doMovie&dontSaveResultsAfterMovie) | (~doMovie)
                if ~exist(saveAddress,'dir')
                    mkdir(saveAddress)
                end
                if ~dontsaveanything
                    switch globalArguments.ferretToAnalyze
                        case {'Munster','Cabecou','Banon','Laguiole'}
                            save([saveAddress saveName],'Y_in','Y_out','M2','shifts2','template2','options_nonrigid','pre_normcorre_cropDimensions','usableTrials','b0_in_initialSize','framesIncluded')
                        case {'Bichonnet','Pigouille','PigouilleOBG'}
                            warning(['I am not saving two variables that might be expected, usableTrials or framesIncluded. Both seem like they might need some rethinking for the new format.']);
                            save([saveAddress saveName],'Y_in','Y_out','M2','shifts2','template2','options_nonrigid','pre_normcorre_cropDimensions','b0_in_initialSize')
                            ';';
                    end
                else
                end
            end
            % /[ferretname]_ProcessedData/[outROI_name]/[inROI_name]/[CCs_removed]
            if ~isempty(variablesToKeep)
                clearvars('-except',variablesToKeep{:})
            end
        end
    end
end
end
