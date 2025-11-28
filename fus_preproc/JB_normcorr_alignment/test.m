filterType = 'keep_tissue';
svdsToTest = 70;
svd_combine_number = 10;
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
    case 'butterworth'
      Band = 60; %%% So as far as I can tell, the formula relating axial speed to cutoff is Cutoff (Hz) = (axial speed (mm/s))/20; This means that I don't want any low-pass applied at all if I can help it. 
      BP = num2str(Band);
      [z,p,k] = butter(6, Band'/(500/2),'low'); %%% 500 is fUS framerate. Was hardcoded by celian so probably won't change. But that's what it is. 
      SOS = zp2sos(z,p,k);
  end
  figure;
  
              for i = 1:svd_combine_number:svdsToTest
                IQ = IQblock;
                IQ = IQ(:,26:80,:);
                IQ_signal = IQ(:,:,1:min(500,size(IQ,3))); %%% seems to only take a maximum of 500 frames... Mine will have 300 or 200 always though. 
                [nz, nx, nt] = size(IQ_signal);
                IQ_signal = reshape(IQ_signal, [nz*nx, nt]);
                cov_matrix = IQ_signal'*IQ_signal;
                [Eig_vect, Eig_val]= eig(cov_matrix);
                Eig_vect=fliplr(Eig_vect);
                Eig_val=rot90(Eig_val,2);
                M_ACP = IQ_signal*Eig_vect;    % on obtient les lambda*u
                skipped_eig_val = i:(i+(svd_combine_number-1));%[1:55];% 190:200];

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
                    Doppler_image = zeros([size(mean(abs(IQF).^2,3)) svdsToTest]);
                end
                Doppler_image(:,:,i) = mean(abs(IQF).^2,3); %%% This is when you get the power doppler signal. "Known to be proportional to blood volume". I wonder if I can still use it for my tissue image, but we'll see. 
               
                imagesc(log(Doppler_image(:,:,i)));colormap('hot');title(num2str(skipped_eig_val));pause(1);
                
              end
                ';';