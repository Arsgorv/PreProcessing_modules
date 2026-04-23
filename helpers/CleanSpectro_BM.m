function [Sc,Th,Epoch]=CleanSpectro_BM(S,f,id,f_thr)

% S=tsd avec time and values of Spectrogram
% [Sc,Th,Epoch]=CleanSpectro(S,f,id)
% id from 1 to 9 (the higher the more removal
% BM add: f_thr for focus only in range of interest (will look for noise > f_thr)

Sp=Data(S);
temp=tsd(Range(S)',mean(Sp(:,f>f_thr)')');

m(1)=percentile(Data(temp),99.9);
m(2)=percentile(Data(temp),99.5);
m(3)=percentile(Data(temp),99);
m(4)=percentile(Data(temp),95);
m(5)=percentile(Data(temp),92);
m(6)=percentile(Data(temp),90);
m(7)=percentile(Data(temp),85);
m(8)=percentile(Data(temp),80);
m(9)=percentile(Data(temp),25);


% 
% m(1)=percentile(Data(temp),0.0099);
% m(2)=percentile(Data(temp),0.09);
% m(3)=percentile(Data(temp),0.1);
% m(4)=percentile(Data(temp),1);
% m(5)=percentile(Data(temp),5);
% m(6)=percentile(Data(temp),10);
% m(7)=percentile(Data(temp),15);
% m(8)=percentile(Data(temp),20);



for i=id
    EpochOK{i}=thresholdIntervals(temp,m(i),'Direction','Below');
end
 

try
    id;
    Epoch=EpochOK{id};
    Sc=Restrict(S,Epoch);
    Th=m(id);
catch
    Epoch=EpochOK{5};
    Sc=Restrict(S,Epoch);
    Th=m(5);
end
