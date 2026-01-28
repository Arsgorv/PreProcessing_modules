function mfile = pick_latest_mfile(stimdir)
m = dir(fullfile(stimdir,'*.m'));
if isempty(m)
    error('No *.m found in %s', stimdir);
end
[~,ix] = max([m.datenum]);
mfile = fullfile(m(ix).folder, m(ix).name);
end