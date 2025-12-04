function React_Passive_epoch_preprocessing
    % archive: delete
        sync_behaviour_ephys(session_dlc{sess})
%% Cut data_cat into trials
% Cut into trials
cut_in_trials_all_sess

% cut into trials
[rawdatacut, trial_timings] = cut_into_trials_AB(data_cat, sess, f_trigs, b_trigs, plt, exp_info);
%         [rawdatacut, trial_timings] = cut_into_trials_AB(raw_data, sess, f_trigs, b_trigs, plt, exp_info);

data_cut_in_trials(:, :, :, :, sess) = rawdatacut;
n_trials = size(rawdatacut, 3);
save('data_cut_in_trials_PreExp', 'data_cut_in_trials', 'trial_timings');


end