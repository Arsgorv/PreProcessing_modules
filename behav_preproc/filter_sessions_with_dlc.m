function session_dlc = filter_sessions_with_dlc(sessions)
k = 1;
session_dlc = {};
for c = 1:numel(sessions)
    dlc_path = fullfile(sessions{c}, 'video');
    files = dir(fullfile(dlc_path, '*_filtered.csv'));
    if ~isempty(files)
        session_dlc{k} = sessions{c};
        k = k + 1;
    else
        disp([sessions{c} ' - No DLC found']);
    end
end
session_dlc = session_dlc(:);
end