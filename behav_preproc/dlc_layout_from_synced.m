function [layoutName, idx] = dlc_layout_from_synced(D, camName)
% D is numeric matrix from synchronized_DLC_data_*.csv
% Col1 = time, Col2 = frame index, then triplets x/y/likelihood per bodypart.
nCols = size(D,2);
if nCols < 5
    error('%s: synchronized DLC has too few columns (%d).', camName, nCols);
end

nTrip = nCols - 2;
if mod(nTrip,3) ~= 0
    error('%s: columns do not match (time+frame+3*k). nCols=%d', camName, nCols);
end

nBp = nTrip/3;

if nBp == 34
    layoutName = 'FULL34';
    bp = { ...
        'implant_1','implant_2', ...
        'ear_1','ear_2','ear_3', ...
        'pupil_1','pupil_2','pupil_3','pupil_4','pupil_5','pupil_6','pupil_7','pupil_8', ...
        'eye_1','eye_2','eye_3','eye_4','eye_5','eye_6','eye_7','eye_8', ...
        'cheek_1','cheek_2','cheek_3','cheek_4', ...
        'nostril_1','nostril_2','nostril_3','nostril_4', ...
        'nose_1','nose_2', ...
        'jaw_1', ...
        'tongue_1', ...
        'spout_1' ...
        };
elseif nBp == 16
    layoutName = 'SMALL16';
    bp = { ...
        'pupil_1','pupil_2','pupil_3','pupil_4','pupil_5','pupil_6','pupil_7','pupil_8', ...
        'eye_1','eye_2','eye_3','eye_4','eye_5','eye_6','eye_7','eye_8' ...
        };
else
    error('%s: unsupported DLC layout: nBodyparts=%d (nCols=%d). Add a mapping.', camName, nBp, nCols);
end

% Build index groups (bodypart indices 1..nBp, in the order above)
idx = struct();

idx.pupil  = find(startsWith(bp,'pupil_'));
idx.eye    = find(startsWith(bp,'eye_'));

if any(startsWith(bp,'implant_'))
    idx.implant = find(startsWith(bp,'implant_'));
else
    idx.implant = [];
end
if any(startsWith(bp,'ear_'))
    idx.ear = find(startsWith(bp,'ear_'));
else
    idx.ear = [];
end
if any(startsWith(bp,'cheek_'))
    idx.cheek = find(startsWith(bp,'cheek_'));
else
    idx.cheek = [];
end
if any(startsWith(bp,'nostril_'))
    idx.nostril = find(startsWith(bp,'nostril_'));
else
    idx.nostril = [];
end
if any(startsWith(bp,'nose_'))
    idx.nose = find(startsWith(bp,'nose_'));
else
    idx.nose = [];
end
if any(startsWith(bp,'jaw_'))
    idx.jaw = find(startsWith(bp,'jaw_'));
else
    idx.jaw = [];
end
if any(startsWith(bp,'tongue_'))
    idx.tongue = find(startsWith(bp,'tongue_'));
else
    idx.tongue = [];
end
if any(startsWith(bp,'spout_'))
    idx.spout = find(startsWith(bp,'spout_'));
else
    idx.spout = [];
end

% For FACE, enforce implant must exist (otherwise normalization breaks)
if strcmp(camName,'FACE') && isempty(idx.implant)
    error('FACE: expected FULL34 layout with implant points, but file looks like %s.', layoutName);
end

end
