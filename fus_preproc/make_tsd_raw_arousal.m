function [tsd_raw, exp_info] = make_tsd_raw_arousal(datapath, plt)
% Preprocessing Reactivation Session
% Reads fUS .scan files and outputs raw_data matrix
% by Arsenii Goriachenkov, ENS - PSL, LSP, Paris
% 2022-2026
% github.com/arsgorv

% datapath example:
% 'Z:\Arsenii\React_Passive\Processed_data\Edel\20220517_2_n_S'

if nargin < 2 || isempty(plt)
    plt = 0;
end

%% load data
filenames = [datapath filesep 'fUS' filesep '*.scan'];
listings = dir(filenames);

nFiles = numel(listings);
disp(strcat(num2str(nFiles), 'found for this arousal session'));

file_order = 1:nFiles;

data = cell(nFiles,1);
epochDuration = nan(nFiles,1);
exp_info.size = cell(nFiles,1);
exp_info.phasename = cell(nFiles,1);

for count = 1:nFiles
    file_idx = file_order(count);
    file_n = listings(file_idx).name;
    disp(file_n)
    
    nfo = h5info([datapath filesep 'fUS' filesep file_n]);
    tmp = squeeze(h5read(nfo.Filename,'/Data')); % Load raw data
    
    % Edel / Chabichou: 3D (x,y,time) after permutation
    if contains(datapath, 'Edel') || contains(datapath, 'Chabichou')
        if size(tmp,2) ~= 128
            tmp = permute(tmp, [2 1 3]);
        end
        % Kosichka/Ficello: 4D (x,y,time,slice) after permutation
    elseif contains(datapath, 'Kosichka') || contains(datapath, 'Ficello')
        if size(tmp,1) ~= 102
            tmp = permute(tmp, [3 1 4 2]);
        end
    else
        error('make_raw_data:UnknownAnimal', ...
            'Unknown animal type in datapath: %s', datapath);
    end
    data{count} = tmp;
    epochDuration(count) = size(tmp,3);
    exp_info.size{count} = size(tmp);
    exp_info.phasename{count} = file_n;
end

save([datapath filesep 'fUS' filesep 'exp_info.mat'], ...
    'exp_info', '-v7.3');


%% Build raw_data 
% Build raw_data and tsd per slice for Edel / Chabichou (3D: x,y,time)
maxT = max(epochDuration);
   
% Build raw_data and tsd per slice for Kosichka and Ficello (4D: x,y,time,slice)
if contains(datapath, 'Kosichka') || contains(datapath, 'Ficello')
    
    Nx = size(data{1},1);
    Ny = size(data{1},2);
    Nslices = size(data{1},4);
    
    if Nslices ~= 4
        warning('make_raw_data:UnexpectedSlices', ...
            'Kosichka/Ficello: expected 4 slices, found %d', Nslices);
    end
    
    % build and save raw_data + tsd per slice
    slice_letters = 'ABCD';
    
    tail = regexp(datapath, '[^\\\/]+$', 'match', 'once');
    
    for s = 1:Nslices
        % raw_slice: x x y x time x phase(3)
        raw_slice = nan(Nx, Ny, maxT, nFiles, 'single');
        raw_data = [];
        for count = 1:nFiles
            padLength = maxT - epochDuration(count);
            % select current slice: Nx x Ny x T
            tmp = data{count}(:,:,:,s);
            tmp = tmp(:,:,1:epochDuration(count));
            raw_data = cat(3, raw_data, tmp);
        end
                
        Nt_all = size(raw_data,3);
        Fs = 2.5; % frames per second
        tvec = linspace(0, (Nt_all/Fs)*1e4, Nt_all)';
        
        M_perm = permute(raw_data, [3 1 2]);      % T x X x Y
        data2D = reshape(M_perm, Nt_all, Nx*Ny);    % T x (X*Y)
        
        tsd_raw.data = tsd(tvec, data2D);
        tsd_raw.Nx = Nx;
        tsd_raw.Ny = Ny;
        
        % RP_tail_A / B / C / D
        out_name = ['RP_data_' tail '_slice_' slice_letters(s) '.mat'];
        save([datapath filesep 'fUS' filesep out_name], ...
            'tsd_raw', '-v7.3');
    end
end

%% CONTROL: Plot raw data
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

%% make params
% els = strsplit(listings(1+(sess-1) * 3).name, '_');
% 
% params.animal = els{1};
% params.session = els{2};
% params.phase = els{3};
% params.slice = els{4};
% params.pair = els{5};
% disp('saving params');
% 
% save('params', 'params');
end