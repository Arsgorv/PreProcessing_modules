function make_cat_tsd(datapath, plt)
% datapath: e.g. 'Z:\Arsenii\React_Passive\Processed_data\Edel\20220517_2_n_S'
% plt: 0/1 for plotting in frames_alignment_AG

if nargin < 2 || isempty(plt)
    plt = 0;
end

tail = regexp(datapath, '[^\\\/]+$', 'match', 'once');

sliceLetters = 'ABCD';

for s = 1:length(sliceLetters)
    slice_label = sliceLetters(s);
    rpFile = [datapath filesep 'fUS' filesep 'RP_data_' tail '_slice_' slice_label '.mat'];

    if ~exist(rpFile, 'file')
        % No such slice in this session, skip
        continue;
    end

    disp(['Processing ' rpFile])

    % Expect tsd_raw inside
    S = load(rpFile);
    if ~isfield(S, 'tsd_raw')
        warning('File %s does not contain tsd_raw. Skipping.', rpFile);
        continue;
    end
    tsd_raw = S.tsd_raw;

    % tsd_raw: time x (Nx*Ny)
    D  = Data(tsd_raw.data);      % [T x (Nx*Ny)]
    Nx = tsd_raw.Nx;
    Ny = tsd_raw.Ny;
    T  = size(D,1);

    % movie: [x y time]
    movie_raw = reshape(D', Nx, Ny, T);

    % run alignment for this slice
    [movie_cat, shifts2, template2, options_nonrigid] = ...
        frames_alignment_AG(datapath, movie_raw, plt, s);

    % back to tsd
    D_cat = reshape(permute(movie_cat, [3 1 2]), T, Nx*Ny);
    tvec  = Range(tsd_raw.data);

    cat_tsd.data = tsd(tvec, D_cat);
    cat_tsd.Nx = Nx;
    cat_tsd.Ny = Ny;

    % append to RP_<tail>_X.mat
    save(rpFile, 'cat_tsd', 'shifts2', 'template2', 'options_nonrigid', '-append');
end
end