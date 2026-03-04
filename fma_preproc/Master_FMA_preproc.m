function Master_FMA_preproc(sessions, opts)
% Master_FMA_preproc (RA experiment compatible)
%
% Wave_clus output location (NEW):
%   <datapath>/ephys/FMA_waveclus/<Phase>/segments/<SegName>/wave_clus/...
%
% Also saves merged unit times per phase:
%   <datapath>/ephys/FMA_waveclus/<Phase>/FMA_units_<Phase>.mat
%
% And generates sanity figures across high-level epochs:
%   <datapath>/ephys/FMA_waveclus/sanity/FMA_epoch_sanity.png

if nargin < 2, opts = struct(); end
if ischar(sessions) || isstring(sessions), sessions = {char(sessions)}; end

% defaults
if ~isfield(opts,'phases'), opts.phases = {'PreSleep','Conditioning','PostSleep','PostTest'}; end
if ~isfield(opts,'fma_channels'), opts.fma_channels = 1:64; end
if ~isfield(opts,'numChannelsFallback'), opts.numChannelsFallback = 97; end 
if ~isfield(opts,'FsFallback'), opts.FsFallback = 30000; end
if ~isfield(opts,'force_sort'), opts.force_sort = false; end
if ~isfield(opts,'do_sanity_figs'), opts.do_sanity_figs = true; end

% common ref: two arrays by default
if ~isfield(opts,'groups'),    opts.groups = {1:32, 33:64}; end
if ~isfield(opts,'commonRef'), opts.commonRef = {1:32, 33:64}; end

% add preprocessing + analysis helpers
thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
addpath(genpath(fullfile(thisDir,'preprocessing')));
addpath(genpath(fullfile(thisDir,'analysis')));

TsRate = 1e4;

for s = 1:numel(sessions)
    datapath = sessions{s};
    disp('------------------------------------------')
    disp(['[Master_FMA_preproc] ' datapath])

    wcRoot = fullfile(datapath, 'ephys', 'FMA_waveclus');
    if ~exist(wcRoot,'dir'), mkdir(wcRoot); end

    mf = fullfile(datapath,'analysis','run_manifest_RAexp.mat');
    isRAexp = exist(mf,'file') == 2;

    if isRAexp
        Srun = load(mf,'RUN'); RUN = Srun.RUN;

        for p = 1:numel(opts.phases)
            phase = opts.phases{p};
            segIdx = find(strcmp({RUN.segments.phase}, phase));
            if isempty(segIdx)
                disp(['  no segments for phase ' phase]);
                continue
            end

            phaseDir = fullfile(wcRoot, phase);
            segBase  = fullfile(phaseDir, 'segments');
            if ~exist(segBase,'dir'), mkdir(segBase); end

            disp(['  phase ' phase ' : ' num2str(numel(segIdx)) ' segment(s)'])

            % -------- 1) wave_clus per segment (saved under ephys/FMA_waveclus) --------
            for ii = 1:numel(segIdx)
                seg = RUN.segments(segIdx(ii));
                segFolder = fullfile(datapath,'ephys', seg.name);

                streamRoot = oe_find_stream(segFolder, 'Acquisition_Board');
                if isempty(streamRoot)
                    streamRoot = oe_find_stream(segFolder, 'Rhythm_FPGA');
                end
                if isempty(streamRoot)
                    warning('    no acquisition stream found in %s', segFolder);
                    continue
                end

                [Fs, nCh] = oe_read_stream_info(streamRoot, opts.FsFallback, opts.numChannelsFallback);

                segOut = fullfile(segBase, seg.name);
                if ~exist(segOut,'dir'), mkdir(segOut); end

                hasTimes = waveclus_has_times(segOut);

                if opts.force_sort || ~hasTimes
                    disp(['    sorting: ' seg.name])
                    for g = 1:numel(opts.groups)
                        chLst = intersect(opts.fma_channels, opts.groups{g});
                        if isempty(chLst), continue; end
                        refCh = opts.commonRef{min(g,numel(opts.commonRef))};
                        refCh = intersect(refCh, 1:nCh);

                        openEphys2wave_clus(streamRoot, chLst, nCh, refCh, segOut, Fs);
                    end
                else
                    disp(['    skip sort (exists): ' seg.name])
                end
            end

            % -------- 2) merge to absolute OE timeline (per phase) --------
            Units = [];
            uCount = 0;

            for ii = 1:numel(segIdx)
                seg = RUN.segments(segIdx(ii));
                segOut = fullfile(segBase, seg.name);

                if ~exist(fullfile(segOut),'dir'), continue; end
                [spk, meta] = load_wave_clus(segOut);
                if isempty(spk), continue; end

                for k = 1:size(spk,2)
                    vec_ms = spk(:,k);
                    vec_ms = vec_ms(isfinite(vec_ms));
                    if isempty(vec_ms), continue; end

                    t_s  = (double(vec_ms)/1000) + seg.t0_s;
                    t_ts = t_s * TsRate;

                    uCount = uCount + 1;
                    Units(uCount).phase = phase;
                    Units(uCount).segName = seg.name;
                    Units(uCount).channel = meta(k,1);
                    Units(uCount).cluster = meta(k,2);
                    Units(uCount).spike_s = t_s(:);
                    Units(uCount).spike_ts = t_ts(:);
                    Units(uCount).unit_id = sprintf('%s_%s_ch%02d_cl%02d', phase, seg.name, meta(k,1), meta(k,2));
                end
            end

            if ~exist(phaseDir,'dir'), mkdir(phaseDir); end
            save(fullfile(phaseDir, ['FMA_units_' phase '.mat']), 'Units', '-v7.3');
            disp(['  saved merged units: ' phase])
        end

        % -------- 3) sanity figures across epochs (columns) --------
        if opts.do_sanity_figs
            try
                FMA_epoch_sanity_figures(datapath, wcRoot, RUN, opts);
            catch ME
                warning('FMA_epoch_sanity_figures failed: %s', ME.message);
            end
        end

    else
        % Non-RA-exp fallback: sort first ephys folder
        ephysDir = fullfile(datapath,'ephys');
        D = dir(ephysDir);
        D = D([D.isdir]);
        D = D(~ismember({D.name},{'.','..','LFPData','FMA_waveclus'}));
        if isempty(D), warning('No ephys subfolders in %s', ephysDir); continue; end

        segFolder = fullfile(ephysDir, D(1).name);
        streamRoot = oe_find_stream(segFolder,'Acquisition_Board');
        if isempty(streamRoot), streamRoot = oe_find_stream(segFolder,'Rhythm_FPGA'); end
        if isempty(streamRoot), warning('No stream found in %s', segFolder); continue; end

        [Fs,nCh] = oe_read_stream_info(streamRoot, opts.FsFallback, opts.numChannelsFallback);

        phaseDir = fullfile(wcRoot, 'Session');
        segOut   = fullfile(phaseDir, 'segments', D(1).name);
        if ~exist(segOut,'dir'), mkdir(segOut); end

        for g = 1:numel(opts.groups)
            chLst = intersect(opts.fma_channels, opts.groups{g});
            if isempty(chLst), continue; end
            refCh = opts.commonRef{min(g,numel(opts.commonRef))};
            refCh = intersect(refCh, 1:nCh);
            openEphys2wave_clus(streamRoot, chLst, nCh, refCh, segOut, Fs);
        end
    end
end
end

% ---------------- helpers (subfunctions) ----------------

function tf = waveclus_has_times(segOut)
tf = false;
wc = fullfile(segOut,'wave_clus');
if ~exist(wc,'dir'), return; end

if ~isempty(dir(fullfile(wc,'times_*.mat')))
    tf = true;
    return
end

D = dir(wc);
D = D([D.isdir]);
D = D(~ismember({D.name},{'.','..'}));
for i = 1:numel(D)
    if ~isempty(dir(fullfile(wc, D(i).name, 'times_*.mat')))
        tf = true;
        return
    end
end
end

function streamRoot = oe_find_stream(segFolder, key)
% robust to both:
%   segFolder/recording1/continuous/<stream>
%   segFolder/Record Node*/experiment*/recording*/continuous/<stream>
streamRoot = '';

cand = dir(fullfile(segFolder,'recording*','continuous',[key '*']));
for i = 1:numel(cand)
    p = fullfile(cand(i).folder, cand(i).name);
    if exist(fullfile(p,'continuous.dat'),'file')
        streamRoot = p;
        return
    end
end

cand = dir(fullfile(segFolder,'Record Node*','experiment*','recording*','continuous',[key '*']));
for i = 1:numel(cand)
    p = fullfile(cand(i).folder, cand(i).name);
    if exist(fullfile(p,'continuous.dat'),'file')
        streamRoot = p;
        return
    end
end
end

function [Fs, nCh] = oe_read_stream_info(streamRoot, FsFallback, nChFallback)
Fs = FsFallback;
nCh = nChFallback;

recDir = fileparts(fileparts(streamRoot)); % .../recordingX
oebin = fullfile(recDir,'structure.oebin');
if exist(oebin,'file') ~= 2, return; end

txt = '';
try, txt = fileread(oebin); catch, return; end

[~,streamName] = fileparts(streamRoot);
ix = strfind(txt, streamName);
if isempty(ix), ix = 1; else, ix = ix(1); end
win = txt(ix:min(ix+9000, numel(txt)));

tok = regexp(win, 'sample_rate\"?\s*:\s*([0-9]+\.?[0-9]*)', 'tokens','once');
if ~isempty(tok), Fs = str2double(tok{1}); end

tok = regexp(win, 'num_channels\"?\s*:\s*([0-9]+)', 'tokens','once');
if ~isempty(tok), nCh = str2double(tok{1}); end
end