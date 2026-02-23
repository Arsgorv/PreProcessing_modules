function STRF_calculation(directory, onsetFile)
%% COMPUTE SPECTROGRAMS FROM TORC
TORClocation = '/media/nas8/OB_ferret_AG_BM/Shropshire/TORCs/';
spectrograms = [];
TORCfnlst = dir([TORClocation '*448*v501*.wav']);
TORCnb = length(TORCfnlst);
for torcNum = 1:TORCnb
    [stimuli(torcNum).file,Fs] = audioread([TORClocation TORCfnlst(torcNum).name]);
    [result(torcNum).s,result(torcNum).w,result(torcNum).t] = spectrogram(stimuli(torcNum).file,440,[],[],Fs,'yaxis');
end
specDat.spectrograms = permute(reshape(vertcat(result(:).s),  size(result(1).s,1), length(result),[]), [3,1,2]);
specDat.freqs= result(torcNum).w;
specDat.time_vector = result(torcNum).t';
soundDur = length(stimuli(torcNum).file)/Fs;
% Display the spectrogram
for ii = []
    figure;
    imagesc(T, freqs, 10*log10(S_log + eps));  % Convert magnitude to dB
    set(gca, 'YScale', 'log');  % Set y-axis to logarithmic scale
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Spectrogram with Logarithmic Frequency Axis');
    colorbar;
    caxis([-10 10]);  % Color scale from -10 to 10 dB
end

%% BIN SPK
file_stim = dir([directory '/Stimuli/*_TORCs.mat']);
file_spikes = dir([directory '/wave_clus/shropshire*']);
TORCmatfile = [file_stim.folder '/' file_stim.name];
spikeFolder{1} = [file_spikes(1).folder '/' file_spikes(1).name '/'];
spikeFolder{2} = [file_spikes(2).folder '/' file_spikes(2).name '/'];

specDat_t = specDat; spectDat_t.time_vector = (0:0.001:3)';
options.plotMe = 0; options.usefirstcycle = 1; options.onlySNR = 1;

for it = 1:2
    [binnedspk{it},chLST{it},clLST{it}] = ARS_binspikes(spikeFolder{it},specDat,TORCmatfile,onsetFile);

    %% COMPUTE STRF
    strf{it} = STRFfromTORC(binnedspk{it},specDat);
    
    %% COMPUTE SNR
    [binnedspk_1kHz{it}] = ARS_binspikes(spikeFolder{it},spectDat_t,TORCmatfile,onsetFile);
%     clear spectDat_t
    d{it} = permute(binnedspk_1kHz{it},[4 1 3 2]); 
    clear binnedspk_1kHz
    
    temp{it} = STRFfromTORC_baphy(d{it},options);
    % ephys_TORC2STRF
    strf{it}.snr = temp{it}.snr;
end

%% PLOT STRF
clear it
NPanelPerFigure = [24 21];
cm = [[linspace(0,1,50)' linspace(0,1,50)' ones(50,1)] ; [ones(50,1) linspace(1,0,50)' linspace(1,0,50)']];


for it = 1:2
    clusterN{it} = size(strf{it}.weights,3);
    figureN{it} = ceil(clusterN{it}/NPanelPerFigure(it));
    for incre = 0:(figureN{it}-1)
        figure;
        for figNum = 1:NPanelPerFigure(it)
            subplot(5,5,figNum);
            cnum = incre*NPanelPerFigure(it)+figNum;
            clear weights
            weights = strf{it}.weights(:,:,cnum);
            weights = weights-median(weights(:));
            imagesc(flipud(squeeze(weights)),max(abs(weights(:)))*[-1 1]);
            % imagesc(flipud(rank1approx));
            if figNum==1
                clear n_freq logf_flip
                n_freq = size(strf{it}.weights,1);
                logf_flip = flipud(strf{it}.log_freg');
                set(gca, 'YTick', 1:10:n_freq, 'YTickLabel', round(logf_flip(1:10:n_freq)));
                ylabel('Frequency (Hz)');
                set(gca,'XTickLabel',[]);
            elseif figNum==21
                clear n_lag
                n_lag = size(strf{it}.weights,2);
                set(gca, 'XTick', 1:10:n_lag, 'XTickLabel', round(strf{it}.lags_ms(1:10:n_lag)));
                xlabel('Lag (ms)');
                set(gca,'YTickLabel',[]);
            else set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); end
            ht=title([sprintf('E.%d C.%.d SNR %.2f',chLST{it}(cnum),clLST{it}(cnum),strf{it}.snr(cnum))],'FontSize',6);
            set(gca,'FontSize',6);
            set(ht,'Interpreter','none');
        end
        drawnow;
        colormap(cm)
    end
end
end