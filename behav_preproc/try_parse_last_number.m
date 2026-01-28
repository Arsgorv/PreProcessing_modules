function x = try_parse_last_number(s)
% Parse last numeric token from a string; return NaN if failed.
x = NaN;
tok = regexp(s,'([-+]?\d*\.?\d+([eE][-+]?\d+)?)','match');
if isempty(tok), return; end
x = str2double(tok{end});
end