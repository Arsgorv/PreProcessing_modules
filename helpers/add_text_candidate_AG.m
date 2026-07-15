function candidateText = add_text_candidate_AG(candidateText, x)

if ischar(x)
    candidateText{end+1,1} = x;
elseif isstring(x)
    candidateText{end+1,1} = char(x);
elseif iscell(x)
    for i = 1:numel(x)
        candidateText = add_text_candidate_AG(candidateText, x{i});
    end
end

end