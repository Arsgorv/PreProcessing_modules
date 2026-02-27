function project = detect_project_from_path(datapath, project_hint)
project = '';
if ~isempty(project_hint)
    project = upper(project_hint);
end
if isempty(project)
    if contains(datapath,'React_Passive') || contains(datapath,'RP_') ...
            || contains(datapath,'Edel') || contains(datapath,'Chabichou') || contains(datapath,'Kosichka') || contains(datapath,'Ficello')
        project = 'RP';
    elseif contains(datapath,'React_Active') || contains(datapath,'RA_') ...
            || contains(datapath,'Tvorozhok') || contains(datapath,'Mochi')
        project = 'RA';
    elseif contains(datapath,'Tonotopy') || contains(datapath,'T_')
        project = 'Tonotopy';        
    else
        error('Cannot infer project from datapath: %s', datapath);
    end
end
if ~any(strcmp(project,{'RP','RA', 'Tonotopy'}))
    error('project must be RP, RA or Tonotopy');
end
end