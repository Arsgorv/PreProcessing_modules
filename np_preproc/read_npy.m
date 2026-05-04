function [data, dtype, shape] = read_npy(npyPath)
% Minimal NumPy .npy reader. Supports 1-D arrays of standard dtypes:
% int8/16/32/64, uint8/16/32/64, float32/64, bool. Little-endian.
fid = fopen(npyPath,'r','l');
if fid < 0, error('read_npy:OpenFail','Cannot open %s', npyPath); end
cleaner = onCleanup(@() fclose(fid));

magic = fread(fid,6,'uint8=>char')';
if ~strcmp(magic, char([147 'NUMPY']))
    error('read_npy:BadMagic','Not an NPY file: %s', npyPath);
end
ver = fread(fid,2,'uint8');
if ver(1) == 1
    hlen = fread(fid,1,'uint16');
else
    hlen = fread(fid,1,'uint32');
end
hdr = fread(fid, double(hlen), 'uint8=>char')';

descrTok = regexp(hdr, '''descr''\s*:\s*''([^'']+)''', 'tokens','once');
shapeTok = regexp(hdr, '''shape''\s*:\s*\(([^)]*)\)', 'tokens','once');
if isempty(descrTok) || isempty(shapeTok)
    error('read_npy:BadHeader','Cannot parse %s', npyPath);
end
dtype = descrTok{1};

parts = strsplit(strtrim(shapeTok{1}), ',');
parts = parts(~cellfun(@(p) isempty(strtrim(p)), parts));
shape = cellfun(@(p) str2double(strtrim(p)), parts);
if isempty(shape), shape = 1; end

% Strip endianness/order char
if numel(dtype) >= 2 && any(dtype(1) == '<>|=')
    base = dtype(2:end);
else
    base = dtype;
end
tc = base(1);
nB = str2double(base(2:end));
switch tc
    case 'i', mlType = sprintf('int%d',  nB*8);
    case 'u', mlType = sprintf('uint%d', nB*8);
    case 'f'
        if nB == 4, mlType = 'single'; else, mlType = 'double'; end
    case 'b', mlType = 'logical';
    otherwise, error('read_npy:Dtype','Unsupported dtype %s in %s', dtype, npyPath);
end

n = prod(max(shape,1));
data = fread(fid, n, [mlType '=>' mlType]);
end