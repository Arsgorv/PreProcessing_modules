function Master_FMA_preproc(session)
%{
 %% Done
% if session == 'F04'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241204/shropshire_2024-12-04_12-56-50_fm_torcs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 58, 61:64]; CommonRefChannels_1 = [33, 35:36, 38:40, 46:47, 49:56, 58, 61:64];
%     chLst_2 = [65:73, 80, 83, 86:96]; CommonRefChannels_2 = [65:73, 80, 83, 86:96];
%     disp('running Shropshire_20241204')
% elseif session == 'F05'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241205/shropshire_2024-12-05_17-10-40_fm_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 58, 61:64]; CommonRefChannels_1 = [33, 35:36, 38:40, 46:56, 58, 61:64];
%     chLst_2 = [65:73, 80, 83:96]; CommonRefChannels_2 = [65:72, 83:96];
%     disp('running Shropshire_20241205')
% elseif session == 'F06'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241206/shropshire_2024-12-06_13-20-17_fm_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 58, 61:64]; CommonRefChannels_1 = [33, 35:36, 38:40, 46:47, 49:56, 58, 61:64];
%     chLst_2 = [83, 86:94]; CommonRefChannels_2 = [83, 86:94];
%     disp('running Shropshire_20241206')
% elseif session == 'F09'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241209/shropshire_2024-12-09_13-53-26_fm_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 58, 61:64]; CommonRefChannels_1 = [33, 35:36, 38:40, 46:56, 58, 61:64];
%     chLst_2 = [65:72, 80, 83, 86:96]; CommonRefChannels_2 = [65:72, 80, 83, 86:94];
%     disp('running Shropshire_20241209')
% elseif session == 'F10'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241210/shropshire_2024-12-10_17-15-07_fm_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 36, 38:40, 46:56, 61:64]; CommonRefChannels_1 = [33, 36, 38:40, 46:56, 61:64];
%     chLst_2 = [65:70, 72, 80, 83, 86:96]; CommonRefChannels_2 = [65:70, 72, 80, 83, 86:96];
%     disp('running Shropshire_20241210')
% elseif session == 'F11'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241211/shropshire_2024-12-11_15-48-50_fm_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 58, 61:64]; CommonRefChannels_1 = [33, 35:36, 38:40, 46:56, 58, 61:64];
%     chLst_2 = [65:72, 80, 86:96]; CommonRefChannels_2 = [65:72, 80, 86:96];
%     disp('running Shropshire_20241211')
% elseif session == 'F12'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241212/shropshire_2024-12-12_16-01-11_fm_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 61:64]; CommonRefChannels_1 = [33, 35:36, 38:40, 46:56, 61:64];
%     chLst_2 = [65:72, 80, 86:96]; CommonRefChannels_2 = [65:72, 80, 86:96];
%     disp('running Shropshire_20241212')
% elseif session == 'F13'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241213/shropshire_2024-12-13_16-54-42_fm_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 36, 38:40, 46:47, 49:56, 61:64]; CommonRefChannels_1 = [33, 36, 38:40, 46:47, 49:56, 61:64];
%     chLst_2 = [65:72, 80, 86:96]; CommonRefChannels_2 = [65:72, 80, 86:96];
%     disp('running Shropshire_20241213')
% elseif session == 'F14'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241214/shropshire_2024-12-14_14-24-00_fm_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38, 40, 46:51, 53:56, 58, 61:64]; CommonRefChannels_1 = [33, 35:36, 38, 40, 46:51, 53:56, 58, 61:64];
%     chLst_2 = [65:72, 80, 86:96]; CommonRefChannels_2 = [65:72, 80, 86:96];
%     disp('running Shropshire_20241214')
% elseif session == 'F21'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241221/shropshire_2024-12-21_14-56-28_fm_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 58, 61:64]; CommonRefChannels_1 = [33, 35:36, 38:40, 46:56, 58, 61:64];
%     chLst_2 = [65:72, 80, 86:96]; CommonRefChannels_2 = [65:72, 80, 86:96];
%     disp('running Shropshire_20241221')
% elseif session == 'F23'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241223/shropshire_2024-12-23_14-22-07_fm_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 58, 61:64]; CommonRefChannels_1 = [33, 35:36, 38:40, 46:56, 58, 61:64];
%     chLst_2 = [65:72, 86, 88:96]; CommonRefChannels_2 = [65:72, 86, 88:96];
%     disp('running Shropshire_20241223')
% elseif session == 'WCH'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241205/shropshire_2024-12-05_17-10-40_fm_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [74:79]; CommonRefChannels_1 = [74:79];
%     disp('running Shropshire_20241205_weird')
% elseif session == 'MC9'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241209/shropshire_2024-12-09_13-53-26_fm_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_2 = [95,96]; CommonRefChannels_2 = [65:72, 80, 83, 86:94];
%     disp('running Shropshire_20241209')
% elseif session == 'H11'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20241211/shropshire_2024-12-11_11-53-11_hf_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 58, 60:64]; CommonRefChannels_1 = [33, 35:36, 38:40, 46:47, 49:56, 58, 60, 62:64];
%     chLst_2 = [65:73, 80, 83, 86:96]; CommonRefChannels_2 = [65:66, 68:73, 80, 86:96];
%     disp('running Shropshire_20241211')
% elseif session == 'H12'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20241212/shropshire_2024-12-12_10-38-15_hf_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 58, 60:64]; CommonRefChannels_1 = [33, 35:36, 38:40, 46:47, 49:56, 58, 60, 62:64];
%     chLst_2 = [65:73, 80, 83, 86:96]; CommonRefChannels_2 = [65:66, 68:73, 80, 86:96];
%     disp('running Shropshire_20241212')
% elseif session == 'H13'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20241213/shropshire_2024-12-13_12-50-47_hf_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 58, 60:64]; CommonRefChannels_1 = [33, 35:36, 38:40, 46:47, 49:56, 58, 60, 62:64];
%     chLst_2 = [65:73, 80, 83, 86:96]; CommonRefChannels_2 = [65:66, 68:73, 80, 86:96];
%     disp('running Shropshire_20241213')
% elseif session == 'H20'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20241220/shropshire_2024-12-20_10-32-03_hf_TORCs_atropine/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 58, 60:64]; CommonRefChannels_1 = [33, 35:36, 38:39, 46:47, 49:56, 58, 61:64];
%     chLst_2 = [65:73, 80, 83, 86:96]; CommonRefChannels_2 = [65:66, 68:71, 73, 80, 86:96];
%     disp('running Shropshire_20241220')
% elseif session == 'H24'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20241224/Shropshire_2024-12-24_12-08-27_hf_TORCs_saline/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 58, 60:64]; CommonRefChannels_1 = [33, 35:36, 38:39, 46:47, 49:56, 58, 61:64];
%     chLst_2 = [65:73, 80, 83, 86:96]; CommonRefChannels_2 = [65:66, 68:71, 73, 80, 86:96];
%     disp('running Shropshire_20241224')
% elseif session == 'H27'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20241227/shropshire_2024-12-27_13-01-34_hf_TORCs_atropine/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 58, 60:64]; CommonRefChannels_1 = [33, 35:36, 38:39, 46:47, 49:56, 58, 61:64];
%     chLst_2 = [65:73, 80, 83, 86:96]; CommonRefChannels_2 = [65:66, 68:71, 73, 80, 86:96];
%     disp('running Shropshire_20241227')
% elseif session == 'H01'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20250101/shropshire_2025-01-01_14-44-46_hf_TORCs_saline/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 58, 60:64]; CommonRefChannels_1 = [33, 35:36, 38:39, 46:47, 49:56, 58, 61:64];
%     chLst_2 = [65:73, 80, 83, 86:96]; CommonRefChannels_2 = [65:66, 68:71, 73, 80, 86:96];
%     disp('running Shropshire_20250101')
% elseif session == 'H05'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20250105/Shropshire_2025-01-05_13-05-48_hf_TORCs_atropine/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 58, 60:64]; CommonRefChannels_1 = [33, 35:36, 38:39, 46:47, 49:56, 58, 61:64];
%     chLst_2 = [65:73, 80, 83, 86:96]; CommonRefChannels_2 = [65:66, 68:71, 73, 80, 86:96];
%     disp('running Shropshire_20250105')    
% elseif session == 'H08'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20250108/2025-01-08_12-16-56/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 58, 60:64]; CommonRefChannels_1 = [33, 35:36, 38:39, 46:47, 49:56, 58, 61:64];
%     chLst_2 = [65:73, 80, 83, 86:96]; CommonRefChannels_2 = [65:66, 68:71, 73, 80, 86:96];
%     disp('running Shropshire_20250108')    
% % end
% elseif session == 'H05'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20241205/shropshire_2024-12-05_13-12-24_hf_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241205')
% elseif session == 'H09'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20241209/shropshire_2024-12-09_09-47-10_hf_TORCS/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241209')
% elseif session == 'H10'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20241210/shropshire_2024-12-10_13-12-10_hf_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241210')
% elseif session == 'H14'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20241214/shropshire_2024-12-14_10-21-52_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241214')
% elseif session == 'H22'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20241222/shropshire_2024-12-22_13-07-54_hf_TORCS_atropine/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241222')    
% elseif session == 'H25'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20241225/shropshire_2024-12-25_14-33-25_hf_TORCs_atropine/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241225')    
% elseif session == 'H26'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20241226/shropshire_2024-12-26_10-06-28_hf_TORCs_saline/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241226')    
% elseif session == 'H28'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20241228/shropshire_2024-12-28_10-35-37_hf_TORCs_saline/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241228')
% elseif session == 'H31'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20241231/shropshire_2024-12-31_12-25-35_hf_TORCs_atropine/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241231')
% elseif session == 'H03'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20250103/shropshire_2025-01-03_10-05-24_hf_TORCs_saline/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20250103')    
% elseif session == 'H41'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20250104/shropshire_2025-01-04_11-48-59_hf_TORCs_saline/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20250104')   
% elseif session == 'F24'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241205/shropshire_2024-12-05_17-10-40_fm_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     numChannels = 113;
%     chLst_1 = [33, 35:36, 38:40, 46:56, 58, 61:64]; CommonRefChannels_1 = [33, 35:36, 38:40, 46:56, 58, 61:64];
%     chLst_2 = [65:73, 80, 83:96]; CommonRefChannels_2 = [65:72, 83:96];
%     disp('running Shropshire_20241205')
% if session == 'L22'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241122/shropshire_2024-11-22_17-07-41_fm_sounds/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241122')
% elseif session == 'L23'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241123/shropshire_2024-11-23_16-54-36_fm_LSP/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241123')
% elseif session == 'L26'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241126/shropshire_2024-11-26_14-34-30_fm_LSP/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241126')    
% elseif session == 'L30'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241130/shropshire_2024-11-30_16-34-15_fm_LSP/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241130') 
%     
% elseif session == 'L03'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241203/shropshire_2024-12-03_14-05-25_fm_LSP/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241203')
% elseif session == 'L17'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241217/shropshire_2024-12-17_13-52-08_LSP_saline/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241217')
% elseif session == 'L18'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241218/shropshire_2024-12-18_14-30-55_fm_LSP/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241218')
% elseif session == 'L24'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241224/shropshire_2024-12-24_16-16-47_fm_LSP_saline/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241224')
% elseif session == 'L82'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241228/shopshire_2024-12-28_14-18-30_fm_LSP_saline/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20241228')
% elseif session == 'L07'
%     root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20250107/shropshire_2025-01-07_12-06-49_fm_LSP_saline/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
%     disp('running Shropshire_20250107')


% Error in openEphys2wave_lus (line 23) - numSamples
if session == 'L28'
    root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241128/shropshire_2024-11-28_13-39-05_fm_LSP/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
    disp('running Shropshire_20241128')
    elseif session == 'H06'
    root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20241206/shropshire_2024-12-06_09-40-33_hf_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
    disp('running Shropshire_20241206')
end
elseif session == 'L30'
    root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20250103/shropshire_2025-01-03_13-43-16_fm_LSP_saline/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
    disp('running Shropshire_20250103')

elseif session == 'H02'
    root = '~/data5/Arsenii/OBG_AG/Shropshire/head-fixed/Shropshire_20250102/shropshire_2025-01-02_12-06-54_hf_TORCs_atropine/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
    disp('running Shropshire_20250102')
end

%%

numChannels = 113;
chLst_1 = [33, 35:36, 38:40, 46:56, 58, 61:64]; CommonRefChannels_1 = [33, 35:36, 38:40, 46:56, 58, 61:64];
chLst_2 = [65:73, 80, 83, 86:96]; CommonRefChannels_2 = [65:73, 80, 83, 86:96];

if session ~= 'L25'
    disp('running the 1st array')
    tic
    openEphys2wave_clus(root,chLst_1,numChannels,CommonRefChannels_1);
    disp(['it took me ' num2str(round(toc/60)) 'm'])
end

disp('running the 2nd array')
tic
openEphys2wave_clus(root,chLst_2,numChannels,CommonRefChannels_2);
disp(['it took me ' num2str(round(toc/60)) 'm'])

%}

tic/60;
root = 'Z:\Arsenii\React_Passive_ephys\Raw_data\Kiri\Kiri_2026-01-08_17-20-02_test\recording1\continuous\Rhythm_FPGA-100.0';
numChannels = 62;
chLst_1 = [1:32]; CommonRefChannels_1 = [1:32];
openEphys2wave_clus(root,chLst_1,numChannels,CommonRefChannels_1);
disp(toc/60)



end