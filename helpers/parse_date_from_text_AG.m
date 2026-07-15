function [t, note] = parse_date_from_text_AG(txt)

t = NaT;
note = '';

tok = regexp(txt, '(\d{4})[-_]?(\d{2})[-_]?(\d{2})', 'tokens', 'once');

if isempty(tok)
    return
end

y = str2double(tok{1});
m = str2double(tok{2});
d = str2double(tok{3});

if isnan(y) || isnan(m) || isnan(d)
    return
end

t = datetime(y, m, d);
note = datestr(t, 'yyyy-mm-dd');

end