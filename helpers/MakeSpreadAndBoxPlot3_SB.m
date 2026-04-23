function [pval , stats_out] = MakeSpreadAndBoxPlot3_SB(A,Cols,X,Legends,varargin)
% Input variables
% A = your data in a cell format A = {Data1,Data2,Data3}
% Cols = The color for each set of data Cols = {[1 0 0],[0 0 1],[0 1 0]}
% If you don't care about colors just put {} and everything will be grey
% X = the position to plot your data X = [1,2,3]
% If you don't care just put []
% Legends = the identity of your datasets for xlabels Legends = {'bla' 'bla' 'bla'}
% If you don't care about colors just put {}
% ShowPoints = 0 for no points = 1 for points
%eventually formating A to cell



%% code from PlotErrorBarN_KJ
% Parse parameter list
for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        error(['Parameter ' num2str(i+2) ' is not a property.']);
    end
    switch(lower(varargin{i}))
        case 'newfig'
            newfig = varargin{i+1};
            if newfig~=0 && newfig ~=1
                error('Incorrect value for property ''newfig''.');
            end
        case 'paired'
            paired = varargin{i+1};
            if paired~=0 && paired ~=1
                error('Incorrect value for property ''paired''.');
            end
        case 'optiontest'
            optiontest = varargin{i+1};
            if ~isstring_FMAToolbox(optiontest, 'ttest' , 'ranksum')
                error('Incorrect value for property ''optiontest''.');
            end
        case 'horizontal'
            horizontal = varargin{i+1};
            if horizontal~=0 && horizontal ~=1
                error('Incorrect value for property ''Horizontal''.');
            end
        case 'showpoints'
            showpoints = varargin{i+1};
            if showpoints~=0 && showpoints ~=1
                error('Incorrect value for property ''showPoints''.');
            end
        case 'barcolors'
            barcolors = varargin{i+1};
        case 'barwidth'
            barwidth = varargin{i+1};
            if ~isvector(barwidth) || length(barwidth)~=1
                error('Incorrect value for property ''barwidth''.');
            end
        case 'x_data'
            x_data = varargin{i+1};
            if ~isvector(x_data)
                error('Incorrect value for property ''x_data''.');
            end
        case 'columntest'
            column_test = varargin{i+1};
            if size(column_test,2)~=2
                error('Incorrect value for property ''ColumnTest''.');
            end
        case 'noncolumntest'
            noncolumn_test = varargin{i+1};
            if size(noncolumn_test,2)~=2
                error('Incorrect value for property ''NonColumnTest''.');
            end
        case 'showsigstar'
            ShowSigstar = varargin{i+1};
            if ~isstring_FMAToolbox(ShowSigstar, 'none' , 'ns', 'sig', 'all')
                error('Incorrect value for property ''ShowSigstar''.');
            end
        case 'y_lim'
            y_lim = varargin{i+1};
            if ~isvector(y_lim) || length(y_lim)~=2
                error('Incorrect value for property ''y_lim''.');
            end
        case 'size_points'
            size_points = varargin{i+1};
        otherwise
            error(['Unknown property ''' num2str(varargin{i}) '''.']);
    end
end

%check if exist and assign default value if not
if ~exist('newfig','var')
    newfig=1;
end
if ~exist('paired','var')
    paired=1;
end
if ~exist('optiontest','var')
    optiontest='ranksum';
end
if ~exist('horizontal','var')
    horizontal=0;
end
if ~exist('showpoints','var')
    showpoints=1;
end
if ~exist('barcolors','var')
    barcolors=[];
end
if ~exist('barwidth','var')
    barwidth=[];
end
if ~exist('x_data','var')
    x_data=[];
end
if ~exist('ShowSigstar','var')
    ShowSigstar='sig';
end
if ~exist('size_points','var')
    size_points=6;
end



B=A; %data not cell
if ~iscell(A)
    tempA=A; A={};
    for i=1:size(tempA,2)
        A{i}=tempA(:,i);
    end
end
if isempty(x_data)
    x = 1:length(A);
else
    x = x_data;
end


if ~strcmpi(ShowSigstar,'none')
    if ~exist('column_test','var')
        column_test = nchoosek(1:length(A),2);
    end
    if exist('noncolumn_test','var')
        [~,indx] = ismember(noncolumn_test, column_test,'rows');
        column_test(indx,:) = [];
    end
    column_test = mat2cell(column_test,ones(1,size(column_test,1)),2);
else
    column_test = cell(0);
end

if strcmpi(ShowSigstar,'sig')
    ShowNS = 0;
    ShowSig = 1;
elseif strcmpi(ShowSigstar,'ns')
    ShowNS = 1;
    ShowSig = 0;
elseif strcmpi(ShowSigstar,'all')
    ShowNS = 1;
    ShowSig = 1;
end



%% Code from Sophie
if isempty(Cols)
    for i = 1:length(A)
        Cols{i} = [0.6 0.6 0.6];
    end
end

if isempty(X)
    for i = 1:length(A)
        X(i) = i*2;
    end
end

if isempty(Legends)
    for i = 1:length(A)
        Legends{i} = num2str(i);
    end
end

ColFact = [1,0.8];
for k = 1 : length(A)
    if sum(isnan(A{k}))<length(A{k})
        %         a=iosr.statistics.boxPlot(X(k),A{k}(:),'boxColor',Cols{k},'lineColor',[0.95 0.95 0.95],'medianColor','k','boxWidth',0.5,'showOutliers',false);
        %         a=iosr.statistics.boxPlot(X(k),A{k}(:),'boxColor',Cols{k},'lineColor',[0.95 0.95 0.95],'medianColor',Cols{k},'boxWidth',0.4,'showOutliers',false);
        a=iosr.statistics.boxPlot(X(k),A{k}(:),'boxColor',Cols{k},'lineColor','none','medianColor',Cols{k},'boxWidth',0.4,'showOutliers',false);
        a.handles.upperWhiskers.Visible='off';a.handles.upperWhiskerTips.Visible='off';
        a.handles.lowerWhiskers.Visible='off';a.handles.lowerWhiskerTips.Visible='off';
        a.handles.medianLines.LineWidth = 10;
        a.handles.medianLines.XData=a.handles.medianLines.XData+[-.1 .1];
        alpha(.3)
        hold on
        
        if showpoints
            %             handlesplot=plotSpread(A{k}(:),'distributionColors','k','xValues',X(k),'spreadWidth',0.8); hold on;
            %             set(handlesplot{1},'MarkerSize',20)
            %             handlesplot=plotSpread(A{k}(:),'distributionColors',Cols{k}*0.4,'xValues',X(k),'spreadWidth',0.8); hold on;
            %             set(handlesplot{1},'MarkerSize',10)
            %             handlesplot=plotSpread(A{k}(:),'distributionColors','k','xValues',X(k),'spreadWidth',0.8); hold on;
            %             set(handlesplot{1},'MarkerSize',10)
            %                         handlesplot=plotSpread(A{k}(:),'distributionColors',[.5 .5 .5],'xValues',X(k),'spreadWidth',0.8); hold on;
            %             set(handlesplot{1},'MarkerSize',10)
            handlesplot=plotSpread(A{k}(:),'distributionColors',Cols{k}*.5,'xValues',X(k),'spreadWidth',.8); hold on;
            set(handlesplot{1},'MarkerSize',size_points)
        end
    end
    MinA(k) = min(A{k});
    MaxA(k) = max(A{k});
end

if paired
    % Preliminary stuff,
    for k= 1: length(A)-1
        N=size(A{k},2); %eventually formating A to cell
        %plot points, indicating value for each individual
        temp=A{k}(and(~isnan(A{k}) , ~isnan(A{k+1})));
        ind = and(~isnan(A{k}) , ~isnan(A{k+1}));
        try
            line([ones(1,sum(ind))*X(k)+.3 ; ones(1,sum(ind))*X(k+1)-.3],[A{k}(ind)' ; A{k+1}(ind)'], 'color',[.7 .7 .7])
        catch
            line([ones(1,sum(ind))*X(k)+.3 ; ones(1,sum(ind))*X(k+1)-.3],[A{k}(ind) ; A{k+1}(ind)], 'color',[.7 .7 .7])
        end
    end
    %if paired data, data are linked together, for each individual, with line
    M = cell2mat(A);
    if size(M,1)==1
        M=cell2mat(A')';
    end
    
    nb_individual=size(M,1); %number of indivudual
    %     for k=1:nb_individual
    %         idx = find(~isnan(M(k,:)));
    %         for i=1:length(idx)-1
    %             try
    %                 line(0.1+[idx(i) idx(i+1)]', [M(k,idx(i)) M(k,idx(i+1))]', 'color',[0.6 0.6 0.6]);
    %             end
    %         end
    %     end
end





xlim([min(X)-1 max(X)+1])
rg = max(MaxA)-min(MinA);
%ylim([min(MinA)-rg*0.1 max(MaxA)+rg*0.1])
if exist('Legends')
    set(gca,'XTick',X,'XTickLabel',Legends)
    box off
    set(gca,'FontSize',10,'Linewidth',1)
    
    
    
    %% STATISTICAL TESTS
    %default column that are tested statistically
    
    thresh_signif=0.05;  % stat param: threshold for statistical significance
    N=length(A); %number of column
    
    if ~isempty(column_test)
        
        % Rank sum Wilcoxon
        if strcmp(optiontest,'ranksum')
            pval=nan(N,N);
            groups = cell(0);
            stats = [];
            for c=1:length(column_test)
                i = column_test{c}(1);
                j = column_test{c}(2);
                if sum(~isnan(A{i}))>2 && sum(~isnan(A{j}))>2
                    if paired
                        idx=find(~isnan(A{i}) & ~isnan(A{j}));
                        [p,h,statstemp]= signrank(A{i}(idx),A{j}(idx));
                        stats_out(i,j) = statstemp.signedrank;
                        stats_out(j,i) = statstemp.signedrank;
                    else
                        [p,h,statstemp]=ranksum(A{i}(~isnan(A{i})),A{j}(~isnan(A{j})));
                        if isfield(statstemp,'zval')
                            stats_out(i,j) = statstemp.zval;
                            stats_out(j,i) = statstemp.zval;
                        else
                            stats_out(i,j) = statstemp.ranksum;
                            stats_out(j,i) = statstemp.ranksum;
                        end
                    end
                    pval(i,j)=p; pval(j,i)=p;
                    
                    
                    if h==1 && ShowSig==1
                        groups{length(groups)+1}=[x(i) x(j)];
                        stats = [stats p];
                    elseif h==0 && ShowNS==1
                        groups{length(groups)+1}=[x(i) x(j)];
                        stats = [stats p];
                    end
                end
            end
            stats(stats>thresh_signif)=nan;
            if ~horizontal
                %             sigstar(groups,stats)
                sigstar(groups,stats)
            else
                sigstarh(groups,stats)
            end
            % T-test
        elseif strcmp(optiontest,'ttest')
            pval=nan(N,N);
            groups = cell(0);
            stats = [];
            % normality
            for i=1:N
                [hi,~,~]=swtest(A{i},0.05);
                hnorm(i)=hi;
            end
            % normal distributions with same variance
            for i=1:N
                for j=i:N
                    [hi,~,~]=vartest2(A{i},A{j});
                    hvar(i,j)=hi;
                end
            end
            % test
            for c=1:length(column_test)
                i = column_test{c}(1);
                j = column_test{c}(2);
                if hnorm(i)+hnorm(j)==0 && hvar(min(i,j),max(i,j))==0
                    if paired
                        [h,p,~,statstemp]=ttest(A{i},A{j});
                    else
                        [h,p,~,statstemp]=ttest2(A{i},A{j});
                    end
                    pval(i,j)=p; pval(j,i)=p;
                    stats_out(i,j) = statstemp.tstat;
                    stats_out(j,i) = statstemp.tstat;
                    
                    if h==1 && ShowSig==1
                        groups{length(groups)+1}=[x(i) x(j)];
                        stats = [stats p];
                    elseif h==0 && ShowNS==1
                        groups{length(groups)+1}=[x(i) x(j)];
                        stats = [stats p];
                    elseif h==0 && ShowNS==0
                        groups{length(groups)+1}=[x(i) x(j)];
                        stats = [stats p];
                    end
                end
            end
            stats(stats>thresh_signif)=nan;
            if or(stats<.05 , ShowNS==1)
                if ~horizontal
                    sigstar(X,stats)
                else
                    sigstarh(X,stats)
                end
            end
        end
    end
    
    %% TERMINATION
    if ~exist('pval','var')
        pval=nan;
    end
    
    
    
    xtickangle(45);
end