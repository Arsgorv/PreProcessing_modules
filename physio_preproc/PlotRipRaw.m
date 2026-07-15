%PlotRipRaw
%
% [M,T] = PlotRipRaw(LFP, events, durations, cleaning, PlotFigure,NewFig)
%
%
%%INPUT
% LFP in tsd
% events : in seconds
% durations in ms
%
%%OUTPUT
% M :   time, mean, std, stdError
% T :   matrix
%
%


function [M,T] = PlotRipRaw(LFP, events, durations, cleaning, PlotFigure, newfig)


%% Inititiation
if exist('durations','var')
    if isdvector(durations,'#1')
        durations = [-durations durations];
    elseif isdvector(durations,'#2','>0')
        durations(1) = -durations(1);
    end
else
    durations = [-50 50]; %50ms
end
durations = durations/1000; %convert in s

if ~exist('cleaning','var')
    cleaning=0;
end
if ~exist('PlotFigure','var')
    PlotFigure=1;
end
if ~exist('newfig','var')
    newfig=1;
end

%params
samplingRate = round(1/median(diff(Range(LFP,'s'))));
nBins = floor(samplingRate*diff(durations)/2)*2+1;

%signals
rg = Range(LFP)/1e4;
LFP_signal = [rg-rg(1) Data(LFP)];

% Changed by Dima Bryzgalov 22/11/2019 (added <-rg(1)>)
%events
if size(events,2)>2
    events_tmp = events(:,2)-rg(1);
else
    events_tmp = events-rg(1);
end


%% Sync
[r,i] = Sync(LFP_signal,events_tmp,'durations',durations);
T = SyncMap(r,i,'durations',durations,'nbins',nBins,'smooth',0);


%% cleaning
if cleaning
    [h,~] = hist(diff(T),100);
    mh=max(h);
    %id=find(mh>mean(mh)+20*std(mh));
    id=find(mh>2*size(T,1)/3);
    % id=find(sum(T)==0);
    idx=find(diff(id)==1);
    id([idx idx+1])=[];
    
    
    if ~isempty(id)
        disp(['Attention problem ! ',num2str(length(id)),' bin(s) at zero'])
        if id(1)>1&id(end)<size(T,2)
            T(:,id)=(T(:,id-1)+T(:,id+1))/2; %!!!! line modified by GL on the 1/05/2015 : new line above :
        elseif id(end)==size(T,2)
            T(:,id(end))=T(:,id(end)-1);
            T(:,id(1:end-1))=(T(:,id(1:end-1)-1)+T(:,id(1:end-1)+1))/2;
        elseif id(1)==1
            T(:,id(1))=T(:,id(1)+1);
            T(:,id(2:end))=(T(:,id(2:end)-1)+T(:,id(2:end)+1))/2;
        end
    end
    %keyboard
end


%% result

%nbin
if size(T,2)>nBins
    nBins=nBins+1;
elseif size(T,2)<nBins
    nBins=nBins-1;
end

%result

try
    M=[((1:nBins)'-ceil(nBins/2))/nBins*diff(durations)' mean(T)' std(T)' stdError(T)'];
catch
    M=[];
    disp('error')
end


%% Plot
try
    if PlotFigure
        
        if newfig
            figure('Color',[1 1 1])
        else
            figure(gcf), hold on,
        end
        
        subplot(1,2,1), hold on
        plot(((1:nBins)'-ceil(nBins/2))/nBins*diff(durations),T,'k');
        plot(((1:nBins)'-ceil(nBins/2))/nBins*diff(durations),mean(T),'r','linewidth',2);
        plot(((1:nBins)'-ceil(nBins/2))/nBins*diff(durations),mean(T)+std(T),'r--');
        plot(((1:nBins)'-ceil(nBins/2))/nBins*diff(durations),mean(T)-std(T),'r--');
        title(['raw data n=' num2str(size(events,1))])
        xlim(durations)
        
        subplot(1,2,2), hold on
        plot(((1:nBins)'-ceil(nBins/2))/nBins*diff(durations), zscore(T')', 'k');
        plot(((1:nBins)'-ceil(nBins/2))/nBins*diff(durations), mean(zscore(T')'), 'r', 'linewidth',2);
        plot(((1:nBins)'-ceil(nBins/2))/nBins*diff(durations), mean(zscore(T')')+std(zscore(T')'), 'r--');
        plot(((1:nBins)'-ceil(nBins/2))/nBins*diff(durations), mean(zscore(T')')-std(zscore(T')'), 'r--');
        title(['zscore n=' num2str(size(T,1))])
        xlim(durations);
    end
end

end

