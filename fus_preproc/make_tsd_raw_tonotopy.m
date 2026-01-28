function [tsd_raw, exp_info] = make_tsd_raw_tonotopy(datapath, plt)
% Preprocessing Tonotopy Session
% Reads one fUS .scan file and outputs raw_data matrix
% by Arsenii Goriachenkov, ENS - PSL, LSP, Paris
% 2022-2026
% github.com/arsgorv

% datapath example:
% 'Z:\Arsenii\React_Passive\Processed_data\Edel\20220517_2_n_S'


if nargin < 2 || isempty(plt)
    plt = 0;
end

%% Load data
filenames = [datapath filesep 'fUS' filesep '*.scan'];
listings = dir(filenames);

if numel(listings) ~= 1
    error('Expected exactly one .scan file in %s', filenames);
end

file_n = listings(1).name;
disp(['Loading: ' file_n])

%% Read data
nfo = h5info([datapath filesep 'fUS' filesep file_n]);
tmp = squeeze(h5read(nfo.Filename,'/Data'));

%% Build raw data
if contains(datapath, 'Edel') || contains(datapath, 'Chabichou')
    % Edel / Chabichou → 3D (x,y,time)
    if size(tmp,2) ~= 128
        tmp = permute(tmp, [2 1 3]);
    end

    Nx = size(tmp,1);
    Ny = size(tmp,2);
    Nt = size(tmp,3);

    raw_data = tmp;   % Nx × Ny × T 

elseif contains(datapath, 'Kosichka') || contains(datapath, 'Ficello')
    % Kosichka/Ficello → 4D (x,y,time,slice)
    if size(tmp,1) ~= 102
        tmp = permute(tmp, [3 1 4 2]);
    end

    Nx = size(tmp,1);
    Ny = size(tmp,2);
    Nt = size(tmp,3);
    Nslices = size(tmp,4);

    if Nslices ~= 4
        warning('Kosichka/Ficello: expected 4 slices, found %d', Nslices);
    end

else
    error('Unknown animal type in datapath: %s', datapath);
end

%% Create tsd_raw (single file)
Fs = 2.5; % Hz
tvec = linspace(0, (Nt/Fs)*1e4, Nt)';

if ndims(tmp) == 3
    % 3D → single slice
    M_perm = permute(raw_data, [3 1 2]);  % T × X × Y
    data2D = reshape(M_perm, Nt, Nx*Ny);

    tsd_raw.data = tsd(tvec, data2D);
    tsd_raw.Nx = Nx;
    tsd_raw.Ny = Ny;

    exp_info.size = size(tmp);
    exp_info.phasename = file_n;

    tail = regexp(datapath, '[^\\\/]+$', 'match', 'once');
    save([datapath filesep 'fUS' filesep 'RP_data_' tail '_single.mat'], ...
         'tsd_raw', '-v7.3');

elseif ndims(tmp) == 4
    % 4D Kosichka → loop over slices A–D
    slice_letters = 'ABCD';
    tail = regexp(datapath, '[^\\\/]+$', 'match', 'once');

    for s = 1:Nslices
        slice_data = tmp(:,:,:,s); % Nx × Ny × T

        M_perm = permute(slice_data, [3 1 2]);  
        data2D = reshape(M_perm, Nt, Nx*Ny);

        tsd_raw.data = tsd(tvec, data2D);
        tsd_raw.Nx = Nx;
        tsd_raw.Ny = Ny;

        out_name = ['RP_data_' tail '_slice_' slice_letters(s) '.mat'];
        save([datapath filesep 'fUS' filesep out_name], ...
             'tsd_raw', '-v7.3');
    end
end

%% Control Plot
if plt
    d = dir([datapath filesep 'fUS' filesep 'RP_*']);
    A = load([datapath filesep 'fUS' filesep d(1).name]);
    dat{1} = Data(A.tsd_raw.data);
    try
        B = load([datapath filesep 'fUS' filesep d(2).name]);
        dat{2} = Data(B.tsd_raw.data);
        
        C = load([datapath filesep 'fUS' filesep d(3).name]);
        dat{3} = Data(C.tsd_raw.data);
        
        D = load([datapath filesep 'fUS' filesep d(4).name]);
        dat{4} = Data(D.tsd_raw.data);
        
    catch
        disp('No B-D slices')
    end
    
    figure; hold on
    
    for k = 1:numel(dat)
        frames = reshape(dat{k}', A.tsd_raw.Nx, A.tsd_raw.Ny, []);
        
        subplot(2, 2, k)
        imagesc(squeeze(mean(frames(:, :, :),3)))
    end
end

end
