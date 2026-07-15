function [tStart, note] = parse_start_datetime_from_expeinfo_AG(datapath)

tStart = NaT;
note = 'start datetime not found in ExpeInfo';

expeFile = fullfile(datapath, 'ExpeInfo.mat');

if ~exist(expeFile, 'file')
    return
end

S = load(expeFile);

if ~isfield(S, 'ExpeInfo')
    note = 'ExpeInfo.mat found but no ExpeInfo variable';
    return
end

ExpeInfo = S.ExpeInfo;

candidateText = {};

if isfield(ExpeInfo, 'PreProcessingInfo')

    P = ExpeInfo.PreProcessingInfo;

    if isfield(P, 'FolderForConcatenation_Ephys')
        candidateText = add_text_candidate_AG(candidateText, P.FolderForConcatenation_Ephys);
    end

    if isfield(P, 'FolderSessionName')
        candidateText = add_text_candidate_AG(candidateText, P.FolderSessionName);
    end

end

if isfield(ExpeInfo, 'SessionType')
    candidateText = add_text_candidate_AG(candidateText, ExpeInfo.SessionType);
end

if isfield(ExpeInfo, 'date')
    candidateText = add_text_candidate_AG(candidateText, ExpeInfo.date);
end

for iText = 1:numel(candidateText)

    thisText = candidateText{iText};

    [thisTime, thisNote] = parse_datetime_from_text_AG(thisText);

    if ~isnat(thisTime)
        tStart = thisTime;
        note = ['start datetime from ExpeInfo: ' thisNote];
        return
    end

end

for iText = 1:numel(candidateText)

    thisText = candidateText{iText};

    [thisDate, thisNote] = parse_date_from_text_AG(thisText);

    if ~isnat(thisDate)
        tStart = thisDate;
        note = ['date only from ExpeInfo: ' thisNote];
        return
    end

end

end