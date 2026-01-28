function cfg = get_trigger_config(datapath)
% get_trigger_config
% Returns a struct describing which LFP channels carry which TTL streams.
%
% Fields (channel indices in LFPData):
%   cfg.fus_ch    : fUS trigger TTL channel (NaN if none)
%   cfg.baphy_ch  : Baphy / trial TTL channel
%   cfg.face_ch   : face camera TTL channel
%   cfg.eye_ch    : eye camera TTL channel
%
% This is shared between React_Passive and React_Active projects.
%
% You should gradually fill in the hard-coded mappings as you stabilise your
% acquisition setups. Unknown channels fall back to an interactive input().

cfg.fus_ch   = NaN;
cfg.baphy_ch = NaN;
cfg.face_ch  = NaN;
cfg.eye_ch   = NaN;
cfg.respi    = NaN;
cfg.heart    = NaN;
cfg.OneBox   = NaN;


% Detect animal from datapath
animals = {'Ficello','Kiri','Edel','Chabichou','Kosichka','Tvorozhok','Shropshire','Brynza','Labneh','Mochi','Brayon'};
animal_name = '';
for i = 1:numel(animals)
    if contains(datapath, animals{i})
        animal_name = animals{i};
        break
    end
end

if isempty(animal_name)
    warning('get_trigger_config:UnknownAnimal', ...
        'Could not detect animal from datapath: %s. Using interactive channel selection.', datapath);
else
    disp(['Detected animal: ' animal_name])
end

% -------------------------------------------------------------------------
% Hard-coded defaults where you are sure of the wiring
% -------------------------------------------------------------------------

switch animal_name
    case {'Edel'}
        cfg.fus_ch   = 37;
        cfg.baphy_ch = 36;
        cfg.face_ch  = [];
        cfg.eye_ch   = [];
        cfg.respi    = [];
        cfg.heart    = [];
        cfg.OneBox   = [];
        cfg.baphy_train_gap_s = 0.10;
    case {'Chabichou'}
        cfg.fus_ch   = [];
        cfg.baphy_ch = [];
        cfg.face_ch  = [];
        cfg.eye_ch   = [];
        cfg.respi    = [];
        cfg.heart    = [];  
        cfg.OneBox   = [];
        cfg.baphy_train_gap_s = 0.10;
    case {'Tvorozhok'} % training
        cfg.fus_ch   = [];
        cfg.respi    = 19;
        cfg.heart    = 20;
        cfg.OneBox   = [];
        cfg.baphy_ch = 22;
        cfg.face_ch  = 23;
        cfg.eye_ch   = 24;
        cfg.baphy_train_gap_s = 0.10;
%     case {'Tvorozhok'} % React Active
%         cfg.fus_ch   = [];
%         cfg.baphy_ch = NaN;
%         cfg.face_ch  = NaN;
%         cfg.eye_ch   = NaN;
%         cfg.respi    = NaN;
%         cfg.heart    = NaN;
%         cfg.OneBox   = NaN;         
    case {'Kosichka'} % React Passive
        cfg.fus_ch   = 21;
        cfg.baphy_ch = 20;
        cfg.face_ch  = 22;
        cfg.eye_ch   = [];
        cfg.respi    = [];
        cfg.heart    = 26;   
        cfg.OneBox   = [];        
        cfg.baphy_train_gap_s = 0.10;
    case {'Mochi'} % training
        cfg.fus_ch   = [];
        cfg.respi    = 19;
        cfg.heart    = 20;
        cfg.OneBox   = [];
        cfg.baphy_ch = 22;
        cfg.face_ch  = 23;
        cfg.eye_ch   = 24;
        cfg.baphy_train_gap_s = 0.10;
%     case {'Mochi'} % React Active
%         cfg.fus_ch   = [];
%         cfg.baphy_ch = NaN;
%         cfg.face_ch  = NaN;
%         cfg.eye_ch   = NaN;
%         cfg.respi    = NaN;
%         cfg.heart    = NaN;       
%         cfg.OneBox   = NaN;         
    case {'Brayon'}
        cfg.fus_ch   = [];
        cfg.respi    = 19;
        cfg.heart    = 20;
        cfg.OneBox   = [];
        cfg.baphy_ch = 22;
        cfg.face_ch  = 23;
        cfg.eye_ch   = 24;
        cfg.baphy_train_gap_s = 0.10;
    case {'Ficello'}
        cfg.fus_ch   = 56;
        cfg.baphy_ch = 55;
        cfg.face_ch  = 57;
        cfg.eye_ch   = [];
        cfg.respi    = [];
        cfg.heart    = 61;        
        cfg.baphy_train_gap_s = 0.10;
    case {'Kiri'}
        cfg.fus_ch   = NaN;
        cfg.baphy_ch = NaN;
        cfg.face_ch  = NaN;
        cfg.eye_ch   = NaN;
        cfg.respi    = NaN;
        cfg.heart    = NaN;        
        cfg.baphy_train_gap_s = 0.10;
    otherwise
        cfg.fus_ch   = NaN;
        cfg.baphy_ch = NaN;
        cfg.face_ch  = NaN;
        cfg.eye_ch   = NaN;
        cfg.respi    = NaN;
        cfg.heart    = NaN;        
        cfg.baphy_train_gap_s = NaN;
end

% -------------------------------------------------------------------------
% Interactive fallback for missing channels
% -------------------------------------------------------------------------

if isnan(cfg.fus_ch)
    ans_ch = input('Enter fUS trigger LFP channel number (or [] if none): ');
    if ~isempty(ans_ch)
        cfg.fus_ch = ans_ch;
    end
end

if isnan(cfg.baphy_ch)
    ans_ch = input('Enter Baphy trigger LFP channel number (or [] if none): ');
    if ~isempty(ans_ch)
        cfg.baphy_ch = ans_ch;
    end
end

if isnan(cfg.face_ch)
    ans_ch = input('Enter FACE camera TTL LFP channel number (or [] if none): ');
    if ~isempty(ans_ch)
        cfg.face_ch = ans_ch;
    end
end

if isnan(cfg.eye_ch)
    ans_ch = input('Enter EYE camera TTL LFP channel number (or [] if none): ');
    if ~isempty(ans_ch)
        cfg.eye_ch = ans_ch;
    end
end

if isnan(cfg.respi)
    ans_ch = input('Enter Respiration TTL LFP channel number (or [] if none): ');
    if ~isempty(ans_ch)
        cfg.respi = ans_ch;
    end
end

if isnan(cfg.heart)
    ans_ch = input('Enter Heart TTL LFP channel number (or [] if none): ');
    if ~isempty(ans_ch)
        cfg.heart = ans_ch;
    end
end

if isnan(cfg.OneBox)
    ans_ch = input('Enter OneBox (NP sync)TTL LFP channel number (or [] if none): ');
    if ~isempty(ans_ch)
        cfg.OneBox = ans_ch;
    end
end

% out_file = fullfile(datapath, 'channel_cfg.mat');
% save(out_file, 'cfg');
% 
% disp(['Saved cfg to: ' out_file])
% disp('----------------------------------')
% 
end
