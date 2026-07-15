


%%
cd('/media/nas7/React_Passive_AG/OBG/Labneh/head-fixed/20230227')


for l=[35 25 26] % EMG, OB, HPC, AuCx
    load([pwd filesep 'LFPData/LFP' num2str(l) '.mat'])
    LFP_ferret{l} = LFP;
end
l = 25; LFP_ferret_Fil{l} = FilterLFP(LFP_ferret{l},[20 100],1025);
l = 26; LFP_ferret_Fil{l} = FilterLFP(LFP_ferret{l},[.1 100],1025);
l = 25; LFP_ferret_Fil2{l} = FilterLFP(LFP_ferret{l},[.1 100],1025);
l = 35; LFP_ferret_Fil{l} = FilterLFP(LFP_ferret{l},[.1 10],1025);

bin = 1;

figure
subplot(121)
i=0;
R = Range(LFP_ferret_Fil{35},'s'); D = Data(LFP_ferret_Fil{35});
plot(R(1:bin:end) , D(1:bin:end)-i*4.5e3-8e3 , 'r' , 'LineWidth' , 2), hold on
i=i+1;
D = Data(LFP_ferret_Fil2{25});
plot(R(1:bin:end) , D(1:bin:end)-i*4.5e3 , 'k'), hold on
i=i+1;
D = Data(LFP_ferret_Fil{26});
plot(R(1:bin:end) , D(1:bin:end)-i*4.5e3 , 'k'), hold on
i=i+1;
D = Data(LFP_ferret_Fil{25});
plot(R(1:bin:end) , D(1:bin:end)-i*4.5e3 , 'k'), hold on
xlim([4100 4110]), ylim([-16e3 1e3]), axis off 


subplot(122)
i=0;
R = Range(LFP_ferret_Fil{35},'s'); D = Data(LFP_ferret_Fil{35});
plot(R(1:bin:end) , runmean(D(1:bin:end) , 50)-i*4.5e3-8e3 , 'r' , 'LineWidth' , 2), hold on
i=i+1;
D = Data(LFP_ferret_Fil2{25});
plot(R(1:bin:end) , D(1:bin:end)-i*4.5e3 , 'k'), hold on
i=i+1;
D = Data(LFP_ferret_Fil{26});
plot(R(1:bin:end) , D(1:bin:end)-i*4.5e3 , 'k'), hold on
i=i+1;
D = Data(LFP_ferret_Fil{25});
plot(R(1:bin:end) , D(1:bin:end)-i*4.5e3 , 'k'), hold on
xlim([6678 6688]), ylim([-16e3 1e3]), axis off 













