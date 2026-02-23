function STRF = STRFfromTORC(binnedspk,specDat)
% [time, stim, cluster]
Y = squeeze(mean(binnedspk,3));
[n_tps, n_stim, n_cluster] = size(Y);

% Going to convert the spectrogram to a dB power scale and the frequency axis to a log frequency axis:
% convert to dB power
% [time, frequency, stimuli]
S = log10(abs(specDat.spectrograms));

% convert to log frequency scale
logf = logspace(log10(200), log10(20000), 100);
logf(end) = specDat.freqs(end);
logf_S = ipermute(interp1(specDat.freqs(2:end), permute(S(:,2:end,:), [2, 1, 3]), logf), [2, 1, 3]);
n_freq = length(logf);
sr = round(1/diff(specDat.time_vector([1,2])));

% We do some reformatting of the spectrogram to group the time and stimulus axis together and to add
% lags/delays, since the STRF is just a weighted sum of lagged regressors/features. In this case the
% regressors/features are the energy of each frequency band.

% spectrogram features
% [time, stim, frequency]
F = permute(logf_S, [1, 3, 2]);

% add lags
% [time, stim, frequency, lag]
for lag = 30
    lags = 0:1:lag;
    Fd = add_delays(F, lags);
    n_lag = length(lags);

    % We now use five-fold cross-validation with ridge regularization to fit the model. Note that we always
    % concatenate the time and stimulus dimension and we concatenate the frequency and lag dimension when
    % fitting the model. The most complicated part of the code is keeping track of these dimensions.
    fprintf(strcat([' fit the STRF \n']))
    % folds for cross validation
    % [time, stim]
    n_folds = 5;
    folds = subdivide(n_stim, n_folds);
    folds = repmat(folds, n_tps, 1);

    % unravel/combine/concatenate time-stim and frequency-lag dimensions
    % Fu: [time/stim, freq/lag]
    % Yu: [time/stim, channel]
%     Fu = unravel(unravel(Fd, [1,2]), [2,3]);
%     Yu = unravel(Y, [1,2]);
    Fu = reshape(Fd, [] , size(Fd,3) , size(Fd,4));
    Fu = reshape(Fu, size(Fu,1) , []);
    Yu = reshape(Y, [] , size(Y,3));

    % get STRF as weights from lagged regression
    % W: [feature(offset,frequency/lag), channel]
    % Wr: [frequency, lag, channel]
    W = regress_weights_from_2way_crossval(Fu, Yu, ...
        'method', 'ridge', 'folds', folds(:), 'std_feats', true, 'demean_feats', true);
    Wr = reshape(W(2:end,:), [n_freq, n_lag, n_cluster]);

    STRF = [];
    STRF.W = W;
    STRF.weights = Wr;
    STRF.rank1approx = [];
    STRF.lags_ms = 1000*lags/sr;
    STRF.log_freg = logf;
    STRF.paramstr = strcat([ 'method-', 'ridge', ...
        '_nfolds-', num2str(n_folds), ...
        '_lagms-' num2str(min(STRF.lags_ms)) '-' num2str(max(STRF.lags_ms)), ...
        '_nlags-' num2str(length(STRF.lags_ms)-1)]);
    % plot example STRF
    for i = 1:n_cluster
        [U,S,V] = svd(Wr(:,:,i), 'econ');
        rank1approx = U(:,1)*S(1,1)*V(:,1)';
        STRF.rank1approx(:,:,i) = rank1approx;
    end
end


%% plot example STRF
if 0
    i = 2; % channel
    figure;
    imagesc(flipud(Wr(:,:,i)));
%     imagesc(flipud(rank1approx));
    set(gca, 'XTick', 1:5:n_lag, 'XTickLabel', 1000*(lags(1:5:n_lag))/sr);
    logf_flip = flipud(logf');
    set(gca, 'YTick', 1:5:n_freq, 'YTickLabel', logf_flip(1:5:n_freq));
    xlabel('Lag (ms)');
    ylabel('Frequency (Hz)');
    title(sprintf('Chan %d STRF', i));
end
%% apply STRF to get model-predicted response // measure (not cross-validated) prediction accuracy
% P = [ones(n_tps*n_stim,1), Fu]*W;
% r = fastcorr(P, Yu);
% figure;
% bar(r);
% xlabel('Channel');
% ylabel('Pearson correlation');

%% measure prediction accuracy using nested cross-validation (slow)
% % P: [time/stim, channel]
% P = regress_predictions_from_3way_crossval(unravel(unravel(Fd, [1,2]), [2,3]), unravel(Y, [1,2]), ...
%     'method', 'ridge', 'test_folds', folds(:), 'train_folds', folds(:), 'std_feats', true, 'demean_feats', true);
% r = fastcorr(P, Yu);
% figure;
% bar(r);
% xlabel('Channel');
% ylabel('Pearson correlation');

