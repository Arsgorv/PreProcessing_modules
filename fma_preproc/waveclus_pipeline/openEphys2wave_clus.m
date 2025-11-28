function openEphys2wave_clus(root,chLst,numChannels,CommonRefChannels)
% root = '~/data5/CST/data-source/mozzarella/mozzarella_005_20241113/passive/mozzarella_CST_passive_2024-11-13_17-02-00_tci-slow/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
% chLst = 1:32; numChannels = 43;

% root = '~/data5/Arsenii/OBG_AG/Shropshire/freely-moving/Shropshire_20241205/shropshire_2024-12-05_17-10-40_fm_TORCs/Record Node 101/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
% chLst = 1:32; numChannels = 113;
% openEphys2wave_clus(root,chLst,numChannels,CommonRefChannels);
% Diplays session name
strLocus2 = strfind(root,'Record Node');
strLocus1 = find(root(1:(strLocus2-2))=='/',1,'last');
folderaddress = ['~/' root(strLocus1+1:strLocus2-2) '_ComRef' num2str(CommonRefChannels(1)) '-' num2str(CommonRefChannels(end)) ];
disp(folderaddress);
mkdir(folderaddress);
cd(folderaddress);

datFile = [root 'continuous.dat'];
% numChannels = 43; % Nombre total de canaux
samplingRate = 30000; % Exemple de fréquence d'échantillonnage, adapte selon ton setup
bytesPerSample = 2; % int16 -> 2 octets

% Création de la mémoire map pour charger les données
dataInfo = dir(datFile); % Récupérer la taille du fichier
numSamples = dataInfo.bytes / (numChannels * bytesPerSample); % Nombre total d'échantillons par canal
m = memmapfile(datFile, ...
    'Format', {'int16', [numChannels, numSamples], 'data'}, ...
    'Writable', false);

if nargin<4; CommonRefChannels = []; end
segmentSize = samplingRate*15; % Nombre d'échantillons à traiter par segment
if ~isempty(CommonRefChannels)
    disp('compute common ref...');
    % Préallocation pour la somme et le comptage
    commonRef = zeros(1, size(m.Data.data, 2)); % Taille en fonction de la longueur totale
    nSegments = ceil(size(m.Data.data, 2) / segmentSize);
    
    % Boucle de traitement par segments
    for iSegment = 1:nSegments
        startIdx = (iSegment-1)*segmentSize + 1;
        endIdx = min(iSegment*segmentSize, size(m.Data.data, 2));
        
        % Charger le segment
        dataSegment = m.Data.data(CommonRefChannels, startIdx:endIdx);
        
        % Insère la moyenne
        commonRef(startIdx:endIdx) = mean(dataSegment, 1);
    end

else
    commonRef = 0;
end

for cnum = chLst
    fprintf('\n Saving raw channel', cnum);
    fprintf(' %d ', cnum);

    % Extraction des données du canal spécifié
    data = m.Data.data(cnum, :);
    data = double(data) - commonRef;

    % Ajout des métadonnées pour wave_clus
    par = struct();
    par.sr = samplingRate; % Fréquence d'échantillonnage
    % Ajouter d'autres paramètres spécifiques à wave_clus si nécessaires
    sr = samplingRate;
    % Sauvegarder les données dans un fichier temporaire pour wave_clus
    rawdataSingleChannel = ['C' num2str(cnum) '.mat'];
    save(rawdataSingleChannel, 'data', 'sr','-v7.3');

    fprintf('\n Processing channel', cnum);
    fprintf(' %d ', cnum);
    rawdataSingleChannel = ['C' num2str(cnum) '.mat'];
    spikesSingleChannel = ['C' num2str(cnum) '_spikes.mat'];
%     WAV_CLUS FUNCTIONS
    Get_spikes(rawdataSingleChannel,'parallel',false);
    Do_clustering(spikesSingleChannel,'parallel',false);
end

