function Dir = MergePathForExperiment(Dir1,Dir2)

% example:
% Dir1=PathForExperimentsML('DPCPX');
% Dir2=PathForExperimentsML('LPS');
% Dir = MergePathForExperiment(Dir1,Dir2);
%
% other function
% RestrictPathForExperiment.m IntersectPathForExperiment.m
 
%% evaluate inputs

if nargin<2
    error('Not enough input arguents')
end


%% check Dir1

for i=1:length(Dir1.path)
    
    try temp=Dir1.manipe{i}; catch, temp='ManipeDir1';end
    Dir1.manipe{i}=temp;
    
    try temp=Dir1.CorrecAmpli(i); catch, temp=1;end
    Dir1.CorrecAmpli(i)=temp;
end


%% check Dir2

for i=1:length(Dir2.path)
    
    try temp=Dir2.manipe{i}; catch, temp='ManipeDir2';end
    Dir2.manipe{i}=temp;
    
    try temp=Dir2.CorrecAmpli(i); catch, temp=1;end
    Dir2.CorrecAmpli(i)=temp;
    
end


%% concatenate Dir1 and Dir2


Dir.path=[Dir1.path,Dir2.path];
try
Dir.name=[Dir1.name,Dir2.name];
end
try    
    Dir.group=[Dir1.group,Dir2.group];
end
Dir.manipe=[Dir1.manipe,Dir2.manipe];
Dir.CorrecAmpli=[Dir1.CorrecAmpli,Dir2.CorrecAmpli];
try 
    Dir.delay=[Dir1.delay,Dir2.delay];
end
try 
    Dir.date=[Dir1.date,Dir2.date];
end
try 
    Dir.ExpeInfo=[Dir1.ExpeInfo,Dir2.ExpeInfo];
end

end

