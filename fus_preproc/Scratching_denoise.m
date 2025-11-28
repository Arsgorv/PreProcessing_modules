%% Setting up
animal = {'Edel', 'Chabichou'};
animal_choice = 2;

%% Edel
% slots = {'F','P','G', 'W','D', 'Q','H', 'S','L', 'U', 'R','I', 'T','C','N','J','V', 'M', 'O', 'K', 'E'}; % full ds
% slots = {'L','I','C', 'N','J', 'M', 'O', 'K', 'E'}; % not yet analysed

% slots = {'F','P','G', 'W','D', 'Q','H', 'S', 'V', 'U', 'R', 'T', 'M', 'O', 'J', 'K'}; % processed

%% Chabichou
slots = {'R','S','P','Q','N','O','C','F','G','J','K','H','I','L','M','D','E','T','U'}; 
% Pn is missing. Smth is weird there
% Om is smth weird.
% Start with C

comp = 'icelos'; % 'icelos' or 'biggerguy'

switch comp
    case 'biggerguy'
        prefix = '/mnt/working1/';
        
    case 'icelos'
        prefix = '/home/';
end

datapath = [prefix 'arsenii/data5/Arsenii/React_Passive_AG/fUS/Processed_Data/' animal{animal_choice}];

sessions = {'Morning','Evening'}'; %'Morning' 'Evening'

folder_code = 'Spectral_analysis/scratching_vaso_classification';
figpath = ['/Figs/' folder_code '/'];

% mask_names = {'ACx_mask', 'Hpc_mask', 'tissue_mask', 'out_mask', 'far_out_mask', 'vessel_mask'};
mask_names = {'ACx_mask', 'Hpc_mask', 'tissue_mask', 'far_out_mask', 'vessel_mask'};

% roi_names = {'Cx', 'Hpc', 'Tissue', 'Out', 'Far-out'};
roi_names = {'Cx', 'Hpc', 'Tissue', 'Far-out', 'vessel'};

num_rois = length(roi_names);

%% Set up parameters
% Spectra
movingwin = [120 10]; % window of 120 samples and step size of 10 samples
params.tapers = [3 5];
params.Fs = 2.5;
params.err = 0;
num_rois = 5; % Assuming there are 5 ROIs based on provided data

% Denoise
fband_noise = [0.8 1.2];
smootime = 1;
mergeduration = 25;
dropduration = 10;
strictness_noise = 1.25; % how many sigmas above the mean you take for the threshold
select_noise_roi = 5; % 1 - Cx; 2 - Hpc; 3 - tissue; 4 - Out; 5 - far_out

% Vasomotion
fband_vaso = [0.16 0.22];
strictness_vaso = 1; % how many sigmas above the mean you take for the threshold
select_vaso_roi = 1; % 1 - Cx; 2 - Hpc; 3 - tissue; 4 - Out; 5 - far_out

% Mean spectra and peaks
search_peaks = {[0.1 0.3], [0.13 0.3], [0.08 0.3],[0.08 0.3],[0.08 0.3]};
smoothing_index = 5;

% Plotting
sbplt_count = {1, 13, 25, 3, 15};
plt_colours  = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880];
font_sizes = 8;

%% Do the job
for slot_counter = 18:size(slots, 2)

    %% Load data
    slot_code = slots{slot_counter};
    
    disp(['working on ' slots{slot_counter}])
    addpath(genpath([datapath '/' slots{slot_counter}]))
    cd([datapath '/' slots{slot_counter}]);
    
    parameters = load('params.mat');
    load('exp_info.mat')
    load('data_cat.mat')
    
    %% Check if spectral data already exists and load it if it does:
    try
        load([cd '/SpectralAnalysis.mat'])
        % Raw info
        spectrum = SpectralAnalysis.RawInfo.spectrum;
        coherence = SpectralAnalysis.RawInfo.coherence;
        phase = SpectralAnalysis.RawInfo.phase;
        t = SpectralAnalysis.RawInfo.t;
        f = SpectralAnalysis.RawInfo.f;
        peakFreqs = SpectralAnalysis.RawInfo.peakFreqs;
        
        % Denoise Info
        strictness_noise = SpectralAnalysis.DenoiseInfo.strictness_noise;
        select_noise_roi = SpectralAnalysis.DenoiseInfo.select_noise_roi;
        
        DenoisedEpoch = SpectralAnalysis.DenoiseInfo.DenoisedEpoch;
        
        spectrum_denoised = SpectralAnalysis.DenoiseInfo.spectrum_denoised;
        coherence_denoised = SpectralAnalysis.DenoiseInfo.coherence_denoised;
        phase_denoised = SpectralAnalysis.DenoiseInfo.phase_denoised;
        t_denoised = SpectralAnalysis.DenoiseInfo.t_denoised;
        peakFreqs_denoised = SpectralAnalysis.DenoiseInfo.peakFreqs_denoised;
        
        % Vasomotion Info
        strictness_vaso = SpectralAnalysis.VasoInfo.strictness_vaso;
        select_vaso_roi = SpectralAnalysis.VasoInfo.select_vaso_roi;
        
        denoised_vaso = SpectralAnalysis.VasoInfo.denoised_vaso;
        denoised_no_vaso = SpectralAnalysis.VasoInfo.denoised_no_vaso;
        
        spectrum_vaso = SpectralAnalysis.VasoInfo.spectrum_vaso;
        spectrum_no_vaso = SpectralAnalysis.VasoInfo.spectrum_no_vaso;
        
        coherence_vaso = SpectralAnalysis.VasoInfo.coherence_vaso;
        coherence_no_vaso = SpectralAnalysis.VasoInfo.coherence_no_vaso;
        
        phase_vaso = SpectralAnalysis.VasoInfo.phase_vaso;
        phase_no_vaso = SpectralAnalysis.VasoInfo.phase_no_vaso;
        
        t_vaso = SpectralAnalysis.VasoInfo.t_vaso;
        t_no_vaso = SpectralAnalysis.VasoInfo.t_no_vaso;
        
        peakFreqs_vaso = SpectralAnalysis.VasoInfo.peakFreqs_vaso;
        peakFreqs_no_vaso = SpectralAnalysis.VasoInfo.peakFreqs_no_vaso;
        
        loaded_data_flag = 1;
    catch
        disp('no SpectralAnalysis.mat found, so I am starting a new analysis')
        loaded_data_flag = 0;
    end
    
    %% Plot figures
    if loaded_data_flag == 0
        f1 = figure;
        sgtitle(['Mean power spectra. Denoised. Z-scored: Session ' slot_code], 'FontWeight', 'bold', 'FontSize', 23)
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    end
    
    %% For a given session, calculate the coherence between all ROIs  
    load([pwd '/masks_' parameters.params.slice{1} '_' parameters.params.pair{1}])

    sess_count = 0;
    for sess = 1:2
        clear filtered_data_noise  data_tsd_noise  tEnveloppe_noise  SmoothData_noise  power_data_noise threshold spectrum  t f  freq_range_indices f_range mean_spectra zs_mean_spectra power_range power_range  TotalEpoch  rois_tsd_denoised spectrum_denoised  t_denoised  f_denoised  freq_range_indices_denoised  f_range_denoised  zs_mean_spectra_denoised  mean_spectra_denoised  power_range_denoised  peaks_denoised
        
        %% Plot figures
        f2 = figure;
        sgtitle(['Spectral overview. Denoising: Session ' slot_code '_' num2str(sess)], 'FontWeight', 'bold', 'FontSize', 23)
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        
        f3 = figure;
        sgtitle(['Spectral overview. Vasomotion: Session ' slot_code '_' num2str(sess)], 'FontWeight', 'bold', 'FontSize', 23)
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        
        %% Load session specific mask and data
        if sess == 1
            ACx_mask = eval(['ACx_mask' '_m']);
            Hpc_mask = eval(['Hpc_mask' '_m']);
            far_out_mask = eval(['far_out_mask' '_m']);
            tissue_mask = eval(['tissue_mask' '_m']);
            vessel_mask = eval(['vessel_mask' '_m']);
            if animal_choice == 1
                in_mask = eval(['in_mask' '_m']);
            end
        elseif sess == 2
            ACx_mask = eval(['ACx_mask' '_n']);
            Hpc_mask = eval(['Hpc_mask' '_n']);
            far_out_mask = eval(['far_out_mask' '_n']);
            tissue_mask = eval(['tissue_mask' '_n']);
            vessel_mask = eval(['vessel_mask' '_n']);            
            if animal_choice == 1 
                in_mask = eval(['in_mask' '_n']);
            end
        end
        if animal_choice == 1
            out_mask = ~in_mask;
            data_cat_masks = {ACx_mask, Hpc_mask, tissue_mask, out_mask, far_out_mask};
        else
            data_cat_masks = {ACx_mask, Hpc_mask, tissue_mask, far_out_mask, vessel_mask};
        end
        
        %% Define ROIs
        rois = cell(num_rois, 1);
        
        %% Create filters
        bpFilt_noise = designfilt('bandpassfir', 'FilterOrder', 20, 'CutoffFrequency1', fband_noise(1), 'CutoffFrequency2', fband_noise(2), 'SampleRate', params.Fs);
        bpFilt_vaso = designfilt('bandpassfir', 'FilterOrder', 20, 'CutoffFrequency1', fband_vaso(1), 'CutoffFrequency2', fband_vaso(2), 'SampleRate', params.Fs);
        
        %% Do the raw calculations
        for idx = 1:num_rois
            %% Prepare data
            data_temp = data_cat(:, :, :, sess) .* data_cat_masks{idx};
            rois{idx} = squeeze(mean(data_temp, 1:2)); % Mean over the first two dimensions
            rois_tsd{idx} = tsd(linspace(1, (size(rois{idx}, 1)/params.Fs)*1e4, size(rois{idx}, 1)), rois{idx});
            [coherence{idx}, phase{idx}, ~, spectrum{idx}, ~, t, f] = cohgramc(rois{idx}, rois{idx}, movingwin, params);
            
            sess_time_tp = Range(rois_tsd{idx});
            TotalEpoch = intervalSet(0, sess_time_tp(end));
            
            %% Define noise power
            freq_filt_indices = (f >= fband_noise(1)) & (f <= fband_noise(2));
            f_filt_range = f(freq_filt_indices);
            
            noise_power{idx} = mean(spectrum{idx}(:, freq_filt_indices),2);
            noise_power_tsd{idx} = tsd(t*1e4, noise_power{idx});
            
            %% Define vasomotion power
            freq_vaso_indices = (f >= fband_vaso(1)) & (f <= fband_vaso(2));
            f_vaso_range = f(freq_vaso_indices);
            
            vaso_power{idx} = mean(spectrum{idx}(:, freq_vaso_indices),2);
            vaso_power_tsd{idx} = tsd(t*1e4, vaso_power{idx});
            
            %% Peaks detection fbands
            freq_range_indices = (f >= search_peaks{idx}(1)) & (f <= search_peaks{idx}(2));
            f_range = f(freq_range_indices);
            
            %% Calculate raw spectra
            mean_spectra{idx} = runmean(mean(spectrum{idx}, 1),smoothing_index);
            zs_mean_spectra{idx} = zscore(mean_spectra{idx});
            power_range{idx} = mean_spectra{idx}(freq_range_indices);
            
            %% Calculate raw peaks
            peaks{idx} = max(power_range{idx});
            peakFreqs{sess}{idx} = f_range(find(power_range{idx} == peaks{idx}));
            
        end
        
        %% Check denoise f2
        noise_satisfaction = 'n';
        peaks_sat = 'n';
        while noise_satisfaction == 'n'
            
            try
                strictness_noise = strictness_noise{sess};
                select_noise_roi = select_noise_roi{sess};
            catch
                disp('')
            end
            
            %% Noise episodes calculations
            for idx = 1:num_rois
                %% Not used. Scratching periods defined by raw data power
                %             %Define scratching episodes
                %             filtered_data_noise{idx} = filtfilt(bpFilt_noise, rois{idx});
                %             data_tsd_noise{idx} = tsd(linspace(1, (size(filtered_data_noise{idx}, 1)/params.Fs)*1e4, size(filtered_data_noise{idx}, 1)), filtered_data_noise{idx});
                %             tEnveloppe_noise{idx} = tsd(Range(data_tsd_noise{idx}), abs(hilbert(Data(data_tsd_noise{idx}))) ); %tsd: hilbert transform then enveloppe
                %             SmoothData_noise{idx} = tsd(Range(tEnveloppe_noise{idx}), runmean(Data(tEnveloppe_noise{idx}), ...
                %                 ceil(smootime/median(diff(Range(tEnveloppe_noise{idx},'s'))))));
                %
                %             % power_time_noise = Range(SmoothData_noise, 's');
                %             power_data_noise{idx} = Data(SmoothData_noise{idx});
                %
                %             [Y_noise{idx},X_noise{idx}]=hist(log(power_data_noise{idx}),1000);
                %
                %             threshold{idx} = mean(power_data_noise{idx}) + strictness*std(power_data_noise{idx});
                %
                %             osc_episodes_all_noise{idx} = thresholdIntervals(SmoothData_noise{idx}, threshold{idx}, 'Direction','Above');
                %             osc_episodes_noise{idx} = dropShortIntervals(osc_episodes_all_noise{idx}, dropduration*1e4);
                %             osc_episodes_noise{idx} = mergeCloseIntervals(osc_episodes_noise{idx}, mergeduration*1e4);
                %
                %
                %             osc_episodes_noise_start{idx} = Start(osc_episodes_noise{idx});
                %             osc_episodes_noise_stop{idx} = Stop(osc_episodes_noise{idx});
                
                %% Scratching periods defined by spectra power in fband_noise
                threshold_noise{idx} = mean(noise_power{idx}) + strictness_noise*std(noise_power{idx});
                
                osc_episodes_all_noise{sess}{idx} = thresholdIntervals(noise_power_tsd{idx}, threshold_noise{idx}, 'Direction','Above');
                osc_episodes_noise{sess}{idx} = dropShortIntervals(osc_episodes_all_noise{sess}{idx}, dropduration*1e4);
                osc_episodes_noise{sess}{idx} = mergeCloseIntervals(osc_episodes_noise{sess}{idx}, mergeduration*1e4);
                
                osc_episodes_noise_start{sess}{idx} = Start(osc_episodes_noise{sess}{idx});
                osc_episodes_noise_stop{sess}{idx} = Stop(osc_episodes_noise{sess}{idx});
                
            end
            
            %% Define epochs based on the selected ROI
            % Denoised
            DenoisedEpoch{sess} = TotalEpoch - osc_episodes_noise{sess}{select_noise_roi};
            
            %% Spectra & peaks calculations
            for idx = 1:num_rois
                
                %% Calculate new spectra
                
                % Denoised
                rois_tsd_denoised{idx} = Restrict(rois_tsd{idx}, DenoisedEpoch{sess});
                [coherence_denoised{idx}, phase_denoised{idx},~,spectrum_denoised{idx},~,t_denoised,~] = cohgramc(Data(rois_tsd_denoised{idx}), Data(rois_tsd_denoised{idx}), movingwin, params);
                mean_spectra_denoised{idx} = runmean(mean(spectrum_denoised{idx}, 1),smoothing_index);
                zs_mean_spectra_denoised{idx} = zscore(mean_spectra_denoised{idx});
                power_range_denoised{idx} = mean_spectra_denoised{idx}(freq_range_indices);
                
                %% Calculate denoised peaks
                if peaks_sat == 'n'
                    peaks_denoised{idx} = max(power_range_denoised{idx});
                    
                    peakFreqs_denoised{sess}{idx} = f_range(find(power_range_denoised{idx} == peaks_denoised{idx}));
                    
                end
            end
            
            %% Substitute with loaded data
            if loaded_data_flag == 1
                peakFreqs_denoised = SpectralAnalysis.DenoiseInfo.peakFreqs_denoised;
            end
            
            %% Plot denoising figure
            set(0,'CurrentFigure',f2)
            
            for idx = 1:num_rois
                subplot(9,4,sbplt_count{idx}:(sbplt_count{idx}+1))
                imagesc(t, f, 10*log10(spectrum{idx})'); axis xy; colormap(custom_cmap('redbluedark'))
                caxis([mean(10*log10(spectrum{idx})',1:2) mean(10*log10(spectrum{idx})',1:2)+15])
                
                title(['Spectrogram of ' roi_names{idx}])
                yline([fband_noise(1) fband_noise(2)], 'LineWidth', 1)
                ylim([0.05 1.25])
                xlim([0 3900])
                for m=1:length(osc_episodes_noise_start{sess}{idx})
                    line([osc_episodes_noise_start{sess}{idx}(m)/1e4 osc_episodes_noise_stop{sess}{idx}(m)/1e4],[1.24 1.24],'color','r','linewidth',3);
                end
                
                subplot(9,4,(sbplt_count{idx}+4):(sbplt_count{idx}+5))
                plot(Range(rois_tsd{idx}, 's'), rois{idx}); xlim([0 3900]); ylimval = ylim;
                title('Raw data')
                for m=1:length(osc_episodes_noise_start{sess}{idx})
                    line([osc_episodes_noise_start{sess}{idx}(m)/1e4 osc_episodes_noise_stop{sess}{idx}(m)/1e4],[ylimval(2) ylimval(2)],'color','r','linewidth',3);
                end
                
                subplot(9,4,(sbplt_count{idx}+8):(sbplt_count{idx}+9));
                plot(Range(noise_power_tsd{idx}, 's'), noise_power{idx});  xlim([0 3900]); ylimval = ylim;
                yline([threshold_noise{idx} threshold_noise{idx}], 'color', 'r', 'LineWidth', 0.5)
                title(['Power in the ' num2str(fband_noise(1)) '-' num2str(fband_noise(2)) 'Hz range'])
                for m=1:length(osc_episodes_noise_start{sess}{idx})
                    line([osc_episodes_noise_start{sess}{idx}(m)/1e4 osc_episodes_noise_stop{sess}{idx}(m)/1e4],[ylimval(2) ylimval(2)],'color','r','linewidth',3);
                end
            end
            
            subplot(9,4,[27 31 35]);
            hold on
            for idx = 1:num_rois
                plot(f, zs_mean_spectra{idx}, 'LineWidth', 2)
                xline([peakFreqs{sess}{idx}], 'color', plt_colours(idx, :), 'LineWidth', 1)
            end
            xlim([0.07 0.5])
            if animal_choice == 1
                legend({['Cx = ' num2str(peakFreqs{sess}{1})], '', ['Hpc = ' num2str(peakFreqs{sess}{2})],  '', ['Tissue = ' num2str(peakFreqs{sess}{3})],  '', ['Out = ' num2str(peakFreqs{sess}{4})],  '', ['Far-out = ' num2str(peakFreqs{sess}{5})], ''})
            else
                legend({['Cx = ' num2str(peakFreqs{sess}{1})], '', ['Hpc = ' num2str(peakFreqs{sess}{2})],  '', ['Tissue = ' num2str(peakFreqs{sess}{3})],  '', ['Far-out = ' num2str(peakFreqs{sess}{4})],  '', ['Vessel = ' num2str(peakFreqs{sess}{5})], ''})
            end
            xlabel('Frequency (Hz)');
            ylabel('Power');
            title(['PSD. Z-scored. Based on ' roi_names{select_noise_roi}])
            set(gca, 'FontSize', font_sizes)
            
            
            subplot(9,4,[28 32 36]);
            hold on
            for idx = 1:num_rois
                plot(f, zs_mean_spectra_denoised{idx}, 'LineWidth', 2)
                xline([peakFreqs_denoised{sess}{idx}], 'color', plt_colours(idx, :), 'LineWidth', 1)
            end
            xlim([0.07 0.5])
            if animal_choice == 1
                legend({['Cx = ' num2str(peakFreqs_denoised{sess}{1})], '', ['Hpc = ' num2str(peakFreqs_denoised{sess}{2})],  '', ['Tissue = ' num2str(peakFreqs_denoised{sess}{3})],  '', ['Out = ' num2str(peakFreqs_denoised{sess}{4})],  '', ['Far-out = ' num2str(peakFreqs_denoised{sess}{5})], ''})
            else
                legend({['Cx = ' num2str(peakFreqs_denoised{sess}{1})], '', ['Hpc = ' num2str(peakFreqs_denoised{sess}{2})],  '', ['Tissue = ' num2str(peakFreqs_denoised{sess}{3})],  '', ['Far-out = ' num2str(peakFreqs_denoised{sess}{4})],  '', ['Vessel = ' num2str(peakFreqs_denoised{sess}{5})], ''})

            end
            xlabel('Frequency (Hz)');
            ylabel('Power');
            title(['Denoised PSD. Z-scored. Based on ' roi_names{select_noise_roi}])
            set(gca, 'FontSize', font_sizes)
            hold off
           
            %% Check satisfaction
            
            if loaded_data_flag == 1
                noise_satisfaction = 'y';
            else
                disp(['Selected ROI is: ' roi_names{select_noise_roi}])
                noise_roi_sat = input('Would you like to change ROI selection for noise? y/n: ', 's');
                if noise_roi_sat == 'y'
                    if animal_choice == 1
                        
                        select_noise_roi = input('Select a new referrence ROI: 1 - Cx, 2 - Hpc, 3 - Tissue, 4 - Out, 5 - Far_out: ');
                    else
                        select_noise_roi = input('Select a new referrence ROI: 1 - Cx, 2 - Hpc, 3 - Tissue, 4 - Far_out, 5 - Vessel: ');
                        
                    end
                    
                end
                
                thr_noise_sat = input('Are you satisfied with the noise threshold? y/n: ', 's');
                if thr_noise_sat == 'n'
                    disp(['Current strictness value is: ' num2str(strictness_noise)])
                    strictness_noise = input("Let's update the threshold. Set a new stricteness for the threshold: ");
                end
                
                peaks_sat = input('Any peaks you want to refine (only pay attention to the denoised ones)? y/n: ', 's');
                if peaks_sat == 'y'
                    disp(["Right, let's check all the peaks"])
                    Cx_peak_sat = input('Are you satisfied with the Cx peak? y/n: ', 's');
                    if Cx_peak_sat == 'n'
                        [peakFreqs_denoised{sess}{1}, ~] = ginput(1);
                    end
                    Hpc_peak_sat = input('Are you satisfied with the Hpc peak? y/n: ', 's');
                    if Hpc_peak_sat == 'n'
                        [peakFreqs_denoised{sess}{2}, ~] = ginput(1);
                    end
                    Tissue_peak_sat = input('Are you satisfied with the Tissue peak? y/n: ', 's');
                    if Tissue_peak_sat == 'n'
                        [peakFreqs_denoised{sess}{3}, ~] = ginput(1);
                    end
                    if animal_choice == 1
                        Out_peak_sat = input('Are you satisfied with the Out peak? y/n: ', 's');
                        if Out_peak_sat == 'n'
                            [peakFreqs_denoised{sess}{4}, ~] = ginput(1);
                        end
                        Far_peak_sat = input('Are you satisfied with the Far-out peak? y/n: ', 's');
                        if Far_peak_sat == 'n'
                            [peakFreqs_denoised{sess}{5}, ~] = ginput(1);
                        end
                    else
                        Far_peak_sat = input('Are you satisfied with the Far-out peak? y/n: ', 's');
                        if Far_peak_sat == 'n'
                            [peakFreqs_denoised{sess}{5}, ~] = ginput(1);
                        end
                        
                        Vessel_peak_sat = input('Are you satisfied with the Vessel peak? y/n: ', 's');
                        if Vessel_peak_sat == 'n'
                            [peakFreqs_denoised{sess}{4}, ~] = ginput(1);
                        end
                        
                    end
                end
                
                if thr_noise_sat == 'y' && noise_roi_sat == 'n' && peaks_sat == 'n'
                    noise_satisfaction = 'y';
                end
                
            end
            
        end
        
        %% Check vasomotion f3
        vaso_satisfaction = 'n';
        vaso_peaks_sat = 'n';
        novaso_peaks_sat = 'n';
        while vaso_satisfaction == 'n'
            try
                strictness_vaso = strictness_vaso{sess};
                select_vaso_roi = select_vaso_roi{sess};
            catch
                disp('')
            end

            %% Vaso episodes calculations
            for idx = 1:num_rois
                %% Vasomotion periods defined by spectra power in fband_vaso
                
                threshold_vaso{idx} = mean(vaso_power{idx}) + strictness_vaso*std(vaso_power{idx});
                
                osc_episodes_all_vaso{sess}{idx} = thresholdIntervals(vaso_power_tsd{idx}, threshold_vaso{idx}, 'Direction','Above');
                osc_episodes_vaso{sess}{idx} = dropShortIntervals(osc_episodes_all_vaso{sess}{idx}, dropduration*1e4);
                osc_episodes_vaso{sess}{idx} = mergeCloseIntervals(osc_episodes_vaso{sess}{idx}, mergeduration*1e4);
                
                osc_episodes_vaso_start{sess}{idx} = Start(osc_episodes_vaso{sess}{idx});
                osc_episodes_vaso_stop{sess}{idx} = Stop(osc_episodes_vaso{sess}{idx});
                
            end
            
            %% Define epochs based on the selected ROI
            
            % Denoised vasomotion
            denoised_vaso{sess} = osc_episodes_vaso{sess}{select_vaso_roi}-osc_episodes_noise{sess}{select_noise_roi};
            denoised_vaso_start{sess} = Start(denoised_vaso{sess});
            denoised_vaso_stop{sess} = End(denoised_vaso{sess});
            % Denoised no vasomotion
            denoised_no_vaso{sess} = TotalEpoch - denoised_vaso{sess} - osc_episodes_noise{sess}{select_noise_roi};
            denoised_no_vaso_start{sess} = Start(denoised_no_vaso{sess});
            denoised_no_vaso_stop{sess} = End(denoised_no_vaso{sess});
            
            %% Spectra & peaks calculations
            for idx = 1:num_rois
                
                %% Calculate vaso/no-vaso spectra
                % Vaso
                rois_tsd_vaso{sess}{idx} = Restrict(rois_tsd{idx}, denoised_vaso{sess});
                [coherence_vaso{sess}{idx}, phase_vaso{sess}{idx},~,spectrum_vaso{sess}{idx},~,t_vaso{sess},~] = cohgramc(Data(rois_tsd_vaso{sess}{idx}), Data(rois_tsd_vaso{sess}{idx}), movingwin, params);
                mean_spectra_vaso{sess}{idx} = runmean(mean(spectrum_vaso{sess}{idx}, 1),smoothing_index);
                zs_mean_spectra_vaso{sess}{idx} = zscore(mean_spectra_vaso{sess}{idx});
                power_range_vaso{sess}{idx} = mean_spectra_vaso{sess}{idx}(freq_range_indices);
                
                % No vaso
                rois_tsd_no_vaso{sess}{idx} = Restrict(rois_tsd{idx}, denoised_no_vaso{sess});
                [coherence_no_vaso{sess}{idx}, phase_no_vaso{sess}{idx},~,spectrum_no_vaso{sess}{idx},~,t_no_vaso{sess},~] = cohgramc(Data( rois_tsd_no_vaso{sess}{idx}), Data( rois_tsd_no_vaso{sess}{idx}), movingwin, params);
                mean_spectra_no_vaso{sess}{idx} = runmean(mean(spectrum_no_vaso{sess}{idx}, 1),smoothing_index);
                zs_mean_spectra_no_vaso{sess}{idx} = zscore(mean_spectra_no_vaso{sess}{idx});
                power_range_no_vaso{sess}{idx} = mean_spectra_no_vaso{sess}{idx}(freq_range_indices);
                
                %% Calculate vaso peaks
                if vaso_peaks_sat == 'n'
                    peaks_vaso{sess}{idx} = max(power_range_vaso{sess}{idx});
                    peakFreqs_vaso{sess}{idx} = f_range(find(power_range_vaso{sess}{idx} == peaks_vaso{sess}{idx}));
                end
                
                %% Calculate no-vaso peaks
                if novaso_peaks_sat == 'n'
                    peaks_no_vaso{sess}{idx} = max(power_range_no_vaso{sess}{idx});
                    peakFreqs_no_vaso{sess}{idx} = f_range(find(power_range_no_vaso{sess}{idx} == peaks_no_vaso{sess}{idx}));
                end
                
            end
            
            %% Substitute with loaded data
            if loaded_data_flag == 1
                peakFreqs_vaso = SpectralAnalysis.VasoInfo.peakFreqs_vaso;
                peakFreqs_no_vaso = SpectralAnalysis.VasoInfo.peakFreqs_no_vaso;
            end
            
            %% Plot vasomotion figure
            set(0,'CurrentFigure', f3)
            
            for idx = 1:num_rois
                subplot(9,4,sbplt_count{idx}:(sbplt_count{idx}+1))
                imagesc(t, f, 10*log10(spectrum{idx})'); axis xy; colormap(custom_cmap('redbluedark'))
                caxis([mean(10*log10(spectrum{idx})',1:2) mean(10*log10(spectrum{idx})',1:2)+15])

                title(['Spectrogram of ' roi_names{idx}])
                yline([fband_vaso(1) fband_vaso(2)], 'LineWidth', 1)
                ylim([0.05 0.4])
                xlim([0 3900])
                for m=1:length(denoised_vaso_start{sess})
                    line([denoised_vaso_start{sess}(m)/1e4 denoised_vaso_stop{sess}(m)/1e4],[0.39 0.39],'color','b','linewidth',3);
                end
                for m=1:length(osc_episodes_noise_start{sess}{select_noise_roi})
                    line([osc_episodes_noise_start{sess}{select_noise_roi}(m)/1e4 osc_episodes_noise_stop{sess}{select_noise_roi}(m)/1e4],[0.39 0.39],'color','r','linewidth',3);
                end
                
                subplot(9,4,(sbplt_count{idx}+4):(sbplt_count{idx}+5))
                plot(Range(rois_tsd{idx}, 's'), rois{idx}); xlim([0 3900]); ylimval = ylim;
                title('Raw data')
                for m=1:length(denoised_vaso_start{sess})
                    line([denoised_vaso_start{sess}(m)/1e4 denoised_vaso_stop{sess}(m)/1e4],[ylimval(2) ylimval(2)],'color','b','linewidth',3);
                end
                for m=1:length(osc_episodes_noise_start{sess}{select_noise_roi})
                    line([osc_episodes_noise_start{sess}{select_noise_roi}(m)/1e4 osc_episodes_noise_stop{sess}{select_noise_roi}(m)/1e4],[ylimval(2) ylimval(2)],'color','r','linewidth',3);
                end
                
                subplot(9,4,(sbplt_count{idx}+8):(sbplt_count{idx}+9));
                plot(Range(vaso_power_tsd{idx}, 's'), vaso_power{idx});  xlim([0 3900]); ylimval = ylim;
                yline([threshold_vaso{idx} threshold_vaso{idx}], 'color', 'r', 'LineWidth', 0.5)
                title(['Power in the ' num2str(fband_vaso(1)) '-' num2str(fband_vaso(2)) 'Hz range'])
                for m=1:length(denoised_vaso_start{sess})
                    line([denoised_vaso_start{sess}(m)/1e4 denoised_vaso_stop{sess}(m)/1e4],[ylimval(2) ylimval(2)],'color','b','linewidth',3);
                end
                for m=1:length(osc_episodes_noise_start{sess}{select_noise_roi})
                    line([osc_episodes_noise_start{sess}{select_noise_roi}(m)/1e4 osc_episodes_noise_stop{sess}{select_noise_roi}(m)/1e4],[ylimval(2) ylimval(2)],'color','r','linewidth',3);
                end
                
            end
            
            subplot(9,4,[27 31 35]);
            hold on
            for idx = 1:num_rois
                plot(f, zs_mean_spectra_vaso{sess}{idx}, 'LineWidth', 2)
                xline([peakFreqs_vaso{sess}{idx}], 'color', plt_colours(idx, :), 'LineWidth', 1)
            end
            xlim([0.07 0.5])
            legend({['Cx = ' num2str(peakFreqs_vaso{sess}{1})], '', ['Hpc = ' num2str(peakFreqs_vaso{sess}{2})],  '', ['Tissue = ' num2str(peakFreqs_vaso{sess}{3})],  '', ['Out = ' num2str(peakFreqs_vaso{sess}{4})],  '', ['Far-out = ' num2str(peakFreqs_vaso{sess}{5})], ''})
            xlabel('Frequency (Hz)');
            ylabel('Power');
            title(['Denoised vasomotion PSD. Z-scored. Based on ' roi_names{select_vaso_roi}])
            set(gca, 'FontSize', font_sizes)
            
            
            subplot(9,4,[28 32 36]);
            hold on
            for idx = 1:num_rois
                plot(f, zs_mean_spectra_no_vaso{sess}{idx}, 'LineWidth', 2)
                xline([peakFreqs_no_vaso{sess}{idx}], 'color', plt_colours(idx, :), 'LineWidth', 1)
            end
            xlim([0.07 0.5])
            legend({['Cx = ' num2str(peakFreqs_no_vaso{sess}{1})], '', ['Hpc = ' num2str(peakFreqs_no_vaso{sess}{2})],  '', ['Tissue = ' num2str(peakFreqs_no_vaso{sess}{3})],  '', ['Out = ' num2str(peakFreqs_no_vaso{sess}{4})],  '', ['Far-out = ' num2str(peakFreqs_no_vaso{sess}{5})], ''})
            xlabel('Frequency (Hz)');
            ylabel('Power');
            title(['Denoised no-vasomotion PSD. Z-scored. Based on ' roi_names{select_vaso_roi}])
            set(gca, 'FontSize', font_sizes)
            hold off
            
            %% Check satisfaction
            if loaded_data_flag == 1
                vaso_satisfaction = 'y';
            else
                vaso_thr_sat = input('Are you satisfied with the vaso power threshold? y/n: ', 's');
                if vaso_thr_sat == 'n'
                    disp(['Current strictness value is: ' num2str(strictness_vaso)])
                    strictness_vaso = input("Let's update the threshold. Set a new stricteness for the threshold: ");
                end
                
                vaso_peaks_sat = input('Any vaso peaks you want to refine? y/n: ', 's');
                if vaso_peaks_sat == 'y'
                    disp(["Right, let's check all the peaks"])
                    Cx_peak_sat = input('Are you satisfied with the Cx peak? y/n: ', 's');
                    if Cx_peak_sat == 'n'
                        [peakFreqs_vaso{sess}{1}, ~] = ginput(1);
                    end
                    Hpc_peak_sat = input('Are you satisfied with the Hpc peak? y/n: ', 's');
                    if Hpc_peak_sat == 'n'
                        [peakFreqs_vaso{sess}{2}, ~] = ginput(1);
                    end
                    Tissue_peak_sat = input('Are you satisfied with the Tissue peak? y/n: ', 's');
                    if Tissue_peak_sat == 'n'
                        [peakFreqs_vaso{sess}{3}, ~] = ginput(1);
                    end
                    Out_peak_sat = input('Are you satisfied with the Out peak? y/n: ', 's');
                    if Out_peak_sat == 'n'
                        [peakFreqs_vaso{sess}{4}, ~] = ginput(1);
                    end
                    Far_peak_sat = input('Are you satisfied with the Far-out peak? y/n: ', 's');
                    if Far_peak_sat == 'n'
                        [peakFreqs_vaso{sess}{5}, ~] = ginput(1);
                    end
                end
                
                novaso_peaks_sat = input('Any no-vaso peaks you want to refine? y/n: ', 's');
                if novaso_peaks_sat == 'y'
                    disp(["Right, let's check all the peaks"])
                    Cx_peak_sat = input('Are you satisfied with the Cx peak? y/n: ', 's');
                    if Cx_peak_sat == 'n'
                        [peakFreqs_no_vaso{sess}{1}, ~] = ginput(1);
                    end
                    Hpc_peak_sat = input('Are you satisfied with the Hpc peak? y/n: ', 's');
                    if Hpc_peak_sat == 'n'
                        [peakFreqs_no_vaso{sess}{2}, ~] = ginput(1);
                    end
                    Tissue_peak_sat = input('Are you satisfied with the Tissue peak? y/n: ', 's');
                    if Tissue_peak_sat == 'n'
                        [peakFreqs_no_vaso{sess}{3}, ~] = ginput(1);
                    end
                    Out_peak_sat = input('Are you satisfied with the Out peak? y/n: ', 's');
                    if Out_peak_sat == 'n'
                        [peakFreqs_no_vaso{sess}{4}, ~] = ginput(1);
                    end
                    Far_peak_sat = input('Are you satisfied with the Far-out peak? y/n: ', 's');
                    if Far_peak_sat == 'n'
                        [peakFreqs_no_vaso{sess}{5}, ~] = ginput(1);
                    end
                end
                
                if vaso_thr_sat == 'y' && vaso_peaks_sat == 'n' && novaso_peaks_sat == 'n'
                    vaso_satisfaction = 'y';
                end
            end
        end
        
        %% Plot Mean spectra across dataset f1
        if loaded_data_flag == 0
            yval = nan(3, 2);
            set(0,'CurrentFigure',f1)
            sbpt1 = subplot(2, 3, 1+sess_count);
            hold on
            for idx = 1:2
                plot(f, zs_mean_spectra_denoised{idx}, 'LineWidth', 2)
                xline([peakFreqs_denoised{sess}{idx}], 'color', plt_colours(idx, :), 'LineWidth', 1)
            end
            xlim([0.07 0.3])
            yval(1, 1:2) = ylim;
            legend({['Cx = ' num2str(peakFreqs_denoised{sess}{1})], '', ['Hpc = ' num2str(peakFreqs_denoised{sess}{2})]})
            title([slot_code '_' num2str(sess) ': Full denoised z-scored spectra'])
            xlabel('Frequency (Hz)');
            ylabel('Power');
            set(gca, 'FontSize', font_sizes+5)
            hold off
            
            sbpt2 = subplot(2, 3, 2+sess_count);
            hold on
            for idx = 1:2
                plot(f, zs_mean_spectra_vaso{sess}{idx}, 'LineWidth', 2)
                xline([peakFreqs_vaso{sess}{idx}], 'color', plt_colours(idx, :), 'LineWidth', 1)
            end
            xlim([0.07 0.3])
            yval(2, 1:2) = ylim;
            legend({['Cx = ' num2str(peakFreqs_vaso{sess}{1})], '', ['Hpc = ' num2str(peakFreqs_vaso{sess}{2})]})
            title([slot_code '_' num2str(sess) ': denoised z-scored spectra of vasomotion episodes'])
            xlabel('Frequency (Hz)');
            ylabel('Power');
            set(gca, 'FontSize', font_sizes+5)
            hold off
            
            sbpt3 = subplot(2, 3, 3+sess_count);
            hold on
            for idx = 1:2
                plot(f, zs_mean_spectra_no_vaso{sess}{idx}, 'LineWidth', 2)
                xline([peakFreqs_no_vaso{sess}{idx}], 'color', plt_colours(idx, :), 'LineWidth', 1)
            end
            xlim([0.07 0.3])
            yval(3, 1:2) = ylim;
            legend({['Cx = ' num2str(peakFreqs_no_vaso{sess}{1})], '', ['Hpc = ' num2str(peakFreqs_no_vaso{sess}{2})]})
            title([slot_code '_' num2str(sess) ': denoised z-scored spectra of no-vasomotion episodes'])
            xlabel('Frequency (Hz)');
            ylabel('Power');
            set(gca, 'FontSize', font_sizes+5)
            hold off
            
            linkaxes([sbpt1, sbpt2, sbpt3], 'y');
            ylim([min(yval(:, 1)) max(yval(:, 2))])
            
            sess_count = 3;
        end
        
        %% Save figures f2 & f3
        saveas(f2, [datapath figpath 'PSD_denoising_zscored_' slot_code '_' num2str(sess)], 'svg')
        saveas(f2, [datapath figpath 'PSD_denoising_zscored_' slot_code '_' num2str(sess)], 'png')
        
        saveas(f3, [datapath figpath 'svg/' 'PSD_vaso_denoised_zscored_' slot_code '_'  num2str(sess)], 'svg')
        saveas(f3, [datapath figpath 'PSD_vaso_denoised_zscored_'  slot_code '_' num2str(sess)], 'png')
        
        %% Save importand data
        if loaded_data_flag == 0
            
            % Raw info
            SpectralAnalysis.RawInfo.spectrum{sess} = spectrum;
            SpectralAnalysis.RawInfo.coherence{sess} = coherence;
            SpectralAnalysis.RawInfo.phase{sess} = phase;
            SpectralAnalysis.RawInfo.t{sess} = t;
            SpectralAnalysis.RawInfo.f{sess} = f;
            SpectralAnalysis.RawInfo.peakFreqs{sess} = peakFreqs{sess};
            
            % Denoise Info
            SpectralAnalysis.DenoiseInfo.strictness_noise{sess} = strictness_noise;
            SpectralAnalysis.DenoiseInfo.select_noise_roi{sess} = select_noise_roi;
            
            SpectralAnalysis.DenoiseInfo.DenoisedEpoch{sess} = DenoisedEpoch{sess};
            
            SpectralAnalysis.DenoiseInfo.spectrum_denoised{sess} = spectrum_denoised;
            SpectralAnalysis.DenoiseInfo.coherence_denoised{sess} = coherence_denoised;
            SpectralAnalysis.DenoiseInfo.phase_denoised{sess} = phase_denoised;
            SpectralAnalysis.DenoiseInfo.t_denoised{sess} = t_denoised;
            SpectralAnalysis.DenoiseInfo.peakFreqs_denoised{sess} = peakFreqs_denoised{sess};
            
            % Vasomotion Info
            SpectralAnalysis.VasoInfo.strictness_vaso{sess} = strictness_vaso;
            SpectralAnalysis.VasoInfo.select_vaso_roi{sess} = select_vaso_roi;
            
            SpectralAnalysis.VasoInfo.denoised_vaso{sess} = denoised_vaso{sess};
            SpectralAnalysis.VasoInfo.denoised_no_vaso{sess} = denoised_no_vaso{sess};
            
            SpectralAnalysis.VasoInfo.spectrum_vaso{sess} = spectrum_vaso;
            SpectralAnalysis.VasoInfo.spectrum_no_vaso{sess} = spectrum_no_vaso;
            
            SpectralAnalysis.VasoInfo.coherence_vaso{sess} = coherence_vaso;
            SpectralAnalysis.VasoInfo.coherence_no_vaso{sess} = coherence_no_vaso;
            
            SpectralAnalysis.VasoInfo.phase_vaso{sess} = phase_vaso;
            SpectralAnalysis.VasoInfo.phase_no_vaso{sess} = phase_no_vaso;
            
            SpectralAnalysis.VasoInfo.t_vaso{sess} = t_vaso;
            SpectralAnalysis.VasoInfo.t_no_vaso{sess} = t_no_vaso;
            
            SpectralAnalysis.VasoInfo.peakFreqs_vaso{sess} = peakFreqs_vaso{sess};
            SpectralAnalysis.VasoInfo.peakFreqs_no_vaso{sess} = peakFreqs_no_vaso{sess};
        end
        
    end
    
    %% Save f1 
    if loaded_data_flag == 0
        
        set(0,'CurrentFigure',f1)
        saveas(gcf, [datapath figpath 'svg/''Slot_overview_' slot_code], 'svg')
        saveas(gcf, [datapath figpath 'Slot_overview_' slot_code], 'png')
        close
    end
    
    %% Save important data
    if loaded_data_flag == 0
        
        % Spectra Info
        SpectralAnalysis.SpectraInfo.movingwin = movingwin;
        SpectralAnalysis.SpectraInfo.params = params;
        SpectralAnalysis.SpectraInfo.smootime = smootime;
        SpectralAnalysis.SpectraInfo.mergeduration = mergeduration;
        SpectralAnalysis.SpectraInfo.dropduration = dropduration;
        SpectralAnalysis.DenoiseInfo.fband_noise = fband_noise;
        SpectralAnalysis.VasoInfo.fband_vaso = fband_vaso;
        SpectralAnalysis.roi_names = roi_names;
        
        save('SpectralAnalysis', 'SpectralAnalysis');
    end
end








