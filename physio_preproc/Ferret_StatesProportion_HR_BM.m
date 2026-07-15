%% Don't forget to include noise to wake 

FM_Sess_mouse = GetSleepSessions_Drugs_BM;
FM_Sess_ferret = GetSleepSessions_Ferret_BM;

HR_Sess_mouse = HeadRestraintSess_BM;
HR_Sess_ferret = HeadRestraintSess_Ferret_BM;


Species={'Mouse','Ferret'};
States={'Wake','Sleep','NREM','REM'};

for sp=1:4
    if sp==1
        Sleep_Path = FM_Sess_mouse.path(1);
    elseif sp==2
        Sleep_Path = FM_Sess_ferret.path(3);
    elseif sp==3
        Sleep_Path = HR_Sess_mouse(10);
    elseif sp==4
        Sleep_Path = HR_Sess_ferret(1);
    end
    for indiv=1%:length(Sleep_Path)
        for f=1:length(Sleep_Path{indiv})
            
            %                 load([Sleep_Path{indiv}{f} '/SleepScoring_OBGamma.mat'],'Epoch','Sleep','Wake','SWSEpoch','REMEpoch','SmoothGamma')
            load([Sleep_Path{indiv}{f} '/SleepScoring_OBGamma.mat'],'Epoch','Sleep','Wake','SWSEpoch','REMEpoch')
            
            for states=1:4
                if states==1
                    State = Wake;
                elseif states==2
                    State = Sleep;
                elseif states==3
                    State = SWSEpoch;
                elseif states==4
                    State = REMEpoch;
                end
                
%                                     Smooth_GammA{sp}{indiv}{f} = SmoothGamma;
                
                if states<3
                    State_Prop(sp,indiv,f,states) = sum(DurationEpoch(and(State , Epoch)))./sum(DurationEpoch(Epoch));
                    State_Number(sp,indiv,f,states) = length(Start(and(State , Epoch)))./(sum(DurationEpoch(Epoch))/3600e4);
                else
                    State_Prop(sp,indiv,f,states) = sum(DurationEpoch(and(State , Epoch)))./sum(DurationEpoch(Sleep));
                    State_Number(sp,indiv,f,states) = length(Start(and(State , Epoch)))./(sum(DurationEpoch(Sleep))/3600e4);
                end
                State_MeanDur(sp,indiv,f,states) = nanmedian(DurationEpoch(and(State , Epoch)))/1e4;
                
                try
                    if states<4
                        State = dropShortIntervals(State , 60e4);
                    else
                        State = dropShortIntervals(State , 10e4);
                    end
                    clear St, St = Start(and(State , Epoch))/1e4;
                    State_FirstOnset(sp,indiv,f,states) = St(1);
                end
            end
            disp(Sleep_Path{indiv}{f})
        end
        State_Prop(sp,indiv,squeeze(State_Prop(sp,indiv,:,2))<.3,4)=NaN;
    end
end
State_Prop(State_Prop==0)=NaN;
State_MeanDur(State_MeanDur==0)=NaN;
State_Number(State_Number==0)=NaN;
State_FirstOnset(State_FirstOnset==0)=NaN;



%%
Cols = {[1 .2 .2],[.2 .2 1]};
X = 1:2;
Legends = {'FM','HR'};
NoLegends = {'',''};

figure
for states=1:2
    subplot(4,2,states)
    MakeSpreadAndBoxPlot3_SB({squeeze(squeeze(State_Prop(2,1,:,states))) squeeze(squeeze(State_Prop(4,1,:,states)))},Cols,X,NoLegends,'showpoints',1,'paired',0)
    if states==1, ylabel('proportion'), end
    ylim([0 1.1])
    title(States(states))
    
    subplot(4,2,states+2)
    MakeSpreadAndBoxPlot3_SB({squeeze(squeeze(State_MeanDur(2,1,:,states))) squeeze(squeeze(State_MeanDur(4,1,:,states)))},Cols,X,NoLegends,'showpoints',1,'paired',0)
    if states==1, ylabel('mean dur (s)'), end
    ylim([0 80])
    
    subplot(4,2,states+4)
    MakeSpreadAndBoxPlot3_SB({squeeze(squeeze(State_Number(2,1,:,states))) squeeze(squeeze(State_Number(4,1,:,states)))},Cols,X,NoLegends,'showpoints',1,'paired',0)
    if states==1, ylabel('ep # occurence / recording hours)'), elseif states==3, ylabel('ep # occurence / sleep hours)'), end
    ylim([0 120])
    
    if states>1
        subplot(4,2,states+6)
        MakeSpreadAndBoxPlot3_SB({squeeze(squeeze(State_FirstOnset(2,1,:,states))) squeeze(squeeze(State_FirstOnset(4,1,:,states)))},Cols,X,Legends,'showpoints',1,'paired',0)
        if states==2, ylabel('first onset (s)'), end
    end
end


figure, sp=3; indiv=1;
for i=1:10
    h=histogram(log10(Data(Smooth_GammA{sp}{indiv}{i})),'BinLimits',[2.08 3.3],'NumBins',1e3); % 91=nansum(and(1<Spectro{3},Spectro{3}<8))
    HistData(i,:) = h.Values;
    recording_time(i) = max(Range(Smooth_GammA{sp}{indiv}{i}))/60e4;
end
close


a=jet;
figure
subplot(1,3,1:2)
makepretty_BM
for i=1:10
    plot(linspace(2.08 , 3.3 , 1e3) , runmean(HistData(i,:),20)./sum(HistData(i,:)) , 'LineWidth' , 2 , 'Color', a(6*i,:)), hold on
end

figure
load('/media/nas7/React_Passive_AG/OBG/Labneh/freely-moving/20230309_2/SleepScoring_OBGamma.mat', 'SmoothGamma')
SmoothGamma = Restrict(SmoothGamma,intervalSet(0,1e7));
h=histogram(log10(Data(SmoothGamma)),'BinLimits',[2.08 3.3],'NumBins',1e3);
HistData_fm = h.Values;
close

plot(linspace(2.08 , 3.16 , 1e3) , runmean(HistData_fm,20)./sum(HistData_fm) , 'k' , 'LineWidth' , 5)
box off
xlabel('OB gamma power (log scale)'), ylabel('#')
colormap jet, colorbar

for sp=[2 4]
    for indiv=1%:length(Sleep_Path)
        for f=3:length(Sleep_Path{indiv})
            
            [gamma_thresh , mu1 , mu2 , std1 , std2 , AshD(sp,f)] = GetGammaThresh(Data(Smooth_GammA{sp}{indiv}{f}), 1, 1);
            
        end
    end
end
AshD(AshD==0)=NaN;

subplot(133)
MakeSpreadAndBoxPlot3_SB({AshD(2,:) AshD(4,:)},Cols,X,Legends,'showpoints',1,'paired',0)
ylim([0 4.5]), ylabel('Ahsman D values')


%% toolbox
% Overview of evolution
figure, sp=4; indiv=1; states=2;
subplot(131)
plot(recording_time)
title('recording time')
box off
makepretty_BM

subplot(232)
plot(squeeze(squeeze(State_Prop(sp,indiv,:,2))))
title('sleep propotion')
box off
makepretty_BM

subplot(233)
plot(squeeze(squeeze(State_MeanDur(sp,indiv,:,2))))
title('sleep mean duration')
box off
makepretty_BM

subplot(235)
plot(squeeze(squeeze(State_Number(sp,indiv,:,2))))
title('sleep mean duration')
box off
makepretty_BM

subplot(236)
plot(squeeze(squeeze(State_FirstOnset(sp,indiv,:,2))))
title('sleep first onset')
box off
makepretty_BM


figure
for i=1:10
    subplot(2,5,i)
    plot(Range(Smooth_GammA{sp}{indiv}{i})/60e4 , runmean(Data(Smooth_GammA{sp}{indiv}{i}),5e4))
    xlim([0 250]), ylim([0 1200])
end






%% others



figure, n=4;
for states=1:4
    subplot(4,n,states)
    MakeSpreadAndBoxPlot3_SB({squeeze(squeeze(State_Prop(2,1,:,states))) squeeze(squeeze(State_Prop(4,1,:,states)))},Cols,X,NoLegends,'showpoints',1,'paired',0)
    if states==1, ylabel('proportion'), end
    ylim([0 1.1])
    title(States(states))
    
    subplot(4,n,states+n)
    MakeSpreadAndBoxPlot3_SB({squeeze(squeeze(State_MeanDur(2,1,:,states))) squeeze(squeeze(State_MeanDur(4,1,:,states)))},Cols,X,NoLegends,'showpoints',1,'paired',0)
    if states==1, ylabel('mean dur (s)'), end
    ylim([0 80])
    
    subplot(4,n,states+2*n)
    MakeSpreadAndBoxPlot3_SB({squeeze(squeeze(State_Number(2,1,:,states))) squeeze(squeeze(State_Number(4,1,:,states)))},Cols,X,NoLegends,'showpoints',1,'paired',0)
    if states==1, ylabel('ep # occurence / recording hours)'), elseif states==3, ylabel('ep # occurence / sleep hours)'), end
    ylim([0 300])
    
    if states>1
        subplot(4,n,states+3*n)
        MakeSpreadAndBoxPlot3_SB({squeeze(squeeze(State_FirstOnset(2,1,:,states))) squeeze(squeeze(State_FirstOnset(4,1,:,states)))},Cols,X,Legends,'showpoints',1,'paired',0)
        if states==2, ylabel('first onset (s)'), end
    end
end




figure, n=4;
for states=1:4
    subplot(4,n,states)
    MakeSpreadAndBoxPlot3_SB({squeeze(squeeze(State_Prop(1,1,:,states))) squeeze(squeeze(State_Prop(3,1,:,states)))},Cols,X,NoLegends,'showpoints',1,'paired',0)
    if states==1, ylabel('proportion'), end
    ylim([0 1.1])
    title(States(states))
    
    subplot(4,n,states+n)
    MakeSpreadAndBoxPlot3_SB({squeeze(squeeze(State_MeanDur(1,1,:,states))) squeeze(squeeze(State_MeanDur(3,1,:,states)))},Cols,X,NoLegends,'showpoints',1,'paired',0)
    if states==1, ylabel('mean dur (s)'), end
    ylim([0 80])
    
    subplot(4,n,states+2*n)
    MakeSpreadAndBoxPlot3_SB({squeeze(squeeze(State_Number(1,1,:,states))) squeeze(squeeze(State_Number(3,1,:,states)))},Cols,X,NoLegends,'showpoints',1,'paired',0)
    if states==1, ylabel('ep # occurence / recording hours)'), elseif states==3, ylabel('ep # occurence / sleep hours)'), end
    ylim([0 300])
    
    if states>1
        subplot(4,n,states+3*n)
        MakeSpreadAndBoxPlot3_SB({squeeze(squeeze(State_FirstOnset(1,1,:,states))) squeeze(squeeze(State_FirstOnset(3,1,:,states)))},Cols,X,Legends,'showpoints',1,'paired',0)
        if states==2, ylabel('first onset (s)'), end
    end
end





figure, sp=3; indiv=1; states=2;
subplot(131)
plot(recording_time)
title('recording time')
box off
makepretty_BM

subplot(232)
plot(squeeze(squeeze(State_Prop(sp,indiv,:,2))))
title('sleep propotion')
box off
makepretty_BM

subplot(233)
plot(squeeze(squeeze(State_MeanDur(sp,indiv,:,2))))
title('sleep mean duration')
box off
makepretty_BM

subplot(235)
plot(squeeze(squeeze(State_Number(sp,indiv,:,2))))
title('sleep mean duration')
box off
makepretty_BM

subplot(236)
plot(squeeze(squeeze(State_FirstOnset(sp,indiv,:,2))))
title('sleep first onset')
box off
makepretty_BM


figure
for i=1:10
    subplot(2,5,i)
    plot(Range(Smooth_GammA{sp}{indiv}{i})/60e4 , runmean(Data(Smooth_GammA{sp}{indiv}{i}),5e4))
    xlim([0 250]), ylim([0 1200])
end





