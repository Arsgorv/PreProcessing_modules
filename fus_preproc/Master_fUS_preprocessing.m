function React_Passive_fUS_preprocessing(sessions)
%{
This is a master script for fUS preprocessing in React Passive project
%}

%% Read fUS files and concatenate. Save raw_data ; Save Exp_info
plt = 0;
for sess = 1:numel(sessions)
    disp(['Working on ' sessions{sess}])
    make_tsd_raw(sessions{sess});
end

%% Align frames within one session. Save data_cat
plt = 0;
for sess = 1:numel(sessions)
    disp(['Working on ' sessions{sess}])
    make_cat_tsd(sessions{sess}, plt);
end
    
%% Draw masks
for sess = 1:numel(sessions)
    disp(['Working on ' sessions{sess}]) 
    draw_masks(sessions{sess});
end

%% %%%%%%%%%%%%%%%%%%%%%%%% NOT USED %%%%%%%%%%%%%%%%%%%%%%%%
    
%% NOT USED: FIGURE: Overview of the dataset

% slots = {['C'],['D','E','F','G','H','I'], ['J','K','L','M','N','O'],['P','Q','R','S','T','U']}; %Chabichou
slots = {['F','P','G','W','D','Q'],['H','S','L','U','E','R'],['I','T','C','N','J','V'],['M', 'O', 'K']}; % full ds

figpath = '/Figs/2024/Dataset_overview/';

for i = 1:length(slots)
    f1 = figure('Visible', 'off');
    sgtitle(['Overview of slices ' slots{i} ' for Edel (log10)'], 'FontWeight', 'bold', 'FontSize', 23)
    set(f1, 'Units', 'Normalized', 'Position', [0.0146 0 0.5000 0.9124]);

    count = 1;
    
    for j = 1:length(slots{i})
        
        slot_code = slots{i}(j);
        disp(['working on ' slot_code])
        addpath(genpath([datapath '/' slot_code]))
        cd([datapath '/' slot_code]);
        load('data_cat.mat')
        
        for sess = 1:2
            subplot(3,4,count)
%             subplot(1,2,count)
            imagesc(log10(squeeze(mean(data_cat(:, :, :, sess), 3))))
            count = count + 1;
            axis equal tight
            xticks('')
            yticks('')
            title(['Slice ' slot_code ', Sess #' num2str(sess)])
        end
    end
    saveas(gcf, [datapath figpath 'svg/' 'Overview_log10_' slots{i}], 'svg')
    saveas(gcf, [datapath figpath 'Overview_log10_' slots{i}], 'png')
    close

end

%% NOT USED: Scratching_denoise
Scratching_denoise

%% NOT USED: Reposition_block_images
%{     
    % 20230802 AG: I'm not happy with RepositionBlockImages approach and will try to use Jeffrey's run_normcorre_batch.m
    [data_repo, xoffset, yoffset, xpeak, ypeak, m]  = RepositionBlockImages(data_repo,'First', 20, plt);     
    % find outliers in repositionning
    outlierIdx_end = find(data_repo(end,end,1,:)~=0);
    outlierIdx_start = find(data_repo(end,1,1,:)~=0);
    outlierIdx = [outlierIdx_start outlierIdx_end];
    disp(size(outlierIdx,1))
    
    if size(outlierIdx,1) < 100 % 100 is a random number here distinguish between correct case and others (ugly I know)
        fprintf('There is %d outliers frames after repositionning. Correcting them...\n',length(outlierIdx))
        for outlierI = outlierIdx'
            data_repo(:,:,1,outlierI) = data_repo(:,:,1,outlierI-1);
        end
        
        data_repo = data_repo(1:size(raw_data,1),1:size(raw_data,2),:,:);
    end
    data_cat(:,:,:,sess) = squeeze(data_repo);
%}
    
%% NOT USED: Align two sessions of the same slice (letter)
%{
%Draw the mask for the realignment
disp('time to draw the mask');
figure
for sess = 1:size(data_cat,4)
    Anat(:,:,sess) = sqrt(snm(data_cat(:,:,:,sess),3)); % Here we get rid of the time parameter to get the average activity across session
    It = uint8(255 * mat2gray(squeeze(Anat(:,:,sess))));
    imagesc(It)
    title('Draw approximate mask for realignment')
    Mask(:,:,sess) = roipoly;
end

% for sess = 1:size(data_cat, 4)
%     els = strsplit(listings(1+(sess-1) * 3).name, '_');
%     
%     params.animal{sess} = els{1};
%     params.session{sess} = els{2};
%     params.phase{sess} = els{3};
%     params.slice{sess} = els{4};
%     params.pair{sess} = els{5}; 
% end
% 
% save('params', 'params');

Resp = permute(data_cat, [1 2 4 3]); %Resp is the same as data_repo

% Alignes across two sessions (morning/evening)
disp('Aligning across two sessions (morning/night)');
[data_aligned,anatomy_aligned,mask_aligned] = RepositionImages_new(Resp,Anat,Mask,20,plt,exp_info);
save('data_aligned', 'data_aligned', 'anatomy_aligned','mask_aligned');
save(['REA_AL_B_' slice ], 'params', 'NewMask', 'NewAnat', 'NewResp')
%}

%% NOT USED: denoise with CCA
% The function corrects reactivation sessions by removing components correlated with outside window activity (CCA)
% and by realigning timepoints in time
%
% INPUT
%     - datapath: path to data
%     - slice: slice number     (in character format eg. '01')
%     - pair: sound pair        (in character format eg. 'p4')
%     - repet: pair repetition  (1 or 2)
% OUTPUT
%     - data = Cortex voxels [n_vox * 9000]
%     - save(['REA_' slice '_' pair '_' num2str(repet)], 'D_c','out','Acx', 'data')
%     - out and acx = masks
%     - D_c = corrected data, in and out
% EXAMPLE
% 
% PreprocessSliceAlign_AG(datapath, slice, slot_code)


% % Remove noise with CCA & creates individual sessions files in the form REA_slice.mat
% data_CCA = CCAlign(data_cat, params); %takes around 10 minutes

%% NOT USED: segmentation with NN 
% disp('Starting the segmentation');
% working_folder = [prefix 'arsenii/data5/Arsenii/React_Passive_AG/Florian_legacy/NN_Maryland'];
% load([working_folder '/SegmentationNN_maryland.mat'],'netg'); 
% 
% cm_RdBu = cbrewer('div','RdBu',64,'pchip');
% cm_RdBu = cm_RdBu(end:-1:1,:);
% 
% Anat = sqrt(snm(NewResp,[3 4])); %Here we average across time and sessions. Just the anatomy of the day
% 
% It = uint8(255 * mat2gray(Anat));
% Ct = semanticseg(It,netg);
% 
% figure('Position', [440 277 845 521]);
% 
% It2 = uint8(255 * mat2gray(Anat));
% CBt = labeloverlay(It2, Ct, 'Transparency',.65 );
% I = imfuse(It2,CBt,'montage');
% imagesc(I);
% axis equal tight
% set(gca,'XTick',[],'YTick',[])
% 
% % manual adjustments
% done = 0
% while done ==0
%     figure('units', 'normalized', 'position', [0.1 0 .8 0.9],'Color','w')
%     imagesc(CBt); axis equal tight
%     title('Select voxel to add to Acx')
%     voxs = roipoly;
%     close
%     Ct(voxs) = 'cortex';
%     CBt = labeloverlay(It, Ct);
%     
%     figure('units', 'normalized', 'position', [0.1 0 .8 0.9],'Color','w')
%     imagesc(CBt);axis equal tight
%     title('Select voxels to add to out')
%     voxs = roipoly;
%     Ct(voxs) = 'out';
%     CBt = labeloverlay(It, Ct);
%     close
%     
%     figure('units', 'normalized', 'position', [0 0 .9 0.9],'Color','w')
%     subplot(1,2,1)
%     imagesc(It2);axis equal tight
%     subplot(1,2,2)
%     imagesc(CBt);axis equal tight
%     title('Select voxels to add to hippo')
%     voxs = roipoly;
%     Ct(voxs) = 'hippo';
%     CBt = labeloverlay(It, Ct);
%     close
%     
%     figure('units', 'normalized', 'position', [0.1 0 .8 0.9],'Color','w')
%     imagesc(CBt);
%     title('Are you done? (enter 1 (yes) or 0 (no)')
%     pause(2)
%     done = input('0 or 1: ');
% end
% 
% Acx = (Ct == 'cortex');
% out = (Ct == 'out' | Ct == 'surface' | Ct == 'bloodvessel');
% hippo = (Ct == 'hippo' );
% 
% save(['REA_AL_B_' slice ], 'Ct', 'Acx', 'out', 'hippo', '-append')
% 
% disp('All done, you`re great!');

%% NOT USED: Visualize data

% start = 1;
% 
% figure('Position',[440 443 643 355]);
% dB_gain = 60;
% fig_title = 'Second 0';
% 
% sig = data;
% % sig = (data - snm(data,3))./snm(data,3);
% % sig = Din_c;
% 
% max_image = max(max(abs(sig(:,:))));
% resp = zscore(snm(sig,[1 2]));
% % mv = zscore(corrMvt);
% for ii = start:size(sig, 3)
%     
%     subplot(2,2,1:2)
%     colormap ('hot(256)');
%     imagesc(real(20*log10(sig(:,:,ii)./max_image))), caxis ([-dB_gain 0]);
%     
%     axis equal tight
%     
%     if mod(ii,5)==0
%         fig_title=sprintf('Frame %d, Second %d', ii, ii/5);
%     end
%     title(fig_title);
%     
%     
%     %     subplot(2,2,3)
%     % %     colormap('gray(256)')
%     % %     imagesc(corrFullMvt(:,:,ii));
%     %     axis equal tight
%     
%     subplot(2,2,3:4)
%     range=max(1,ii-50):min(ii+120,length(resp));
%     plot(range,resp(range))
%     hold all
%     %plot(range,bytr(range,:));
%     %plot(range,movPCs(range,4))
%     %     plot(range,mv(range))
%     %     plot(range,zsound(range))
% %     vline(ii)
%     hold off
%     ylim([-2 2])
%     xlim([range(1) range(end)])
%     legend({'Resp in','Resp sub','Resp out'})
%     for j=range
%         if mod(j,29)==0 %12s * 400ms
%             xline(j,'HandleVisibility','off')
%         end
%     end
%     legend({'Resp','Mov','Sound'})
%     drawnow,
%     pause(0.001);
% end

end