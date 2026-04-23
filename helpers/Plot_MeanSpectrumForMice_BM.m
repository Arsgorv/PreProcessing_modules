
% function to plot mice mean spectrums quickly
% varargin can be 'color' and 'threshold', cut off in frequency range for noise detection
% input must be Data lines for mice, frequency in columns
% output is graph parameters
% example: h = Plot_MeanSpectrumForMice_BM(Data , 'color' , 'r' , 'threshold' , 13)


function [h , MaxPowerValues , Freq_Max] = Plot_MeanSpectrumForMice_BM(Data , varargin)

RangeLow = linspace(0.1526,20,261);
RangeHigh = linspace(22,98,32);
RangeVHigh = linspace(22,249,94);
RangeLow2 = linspace(1.0681,20,249);
RangeMiddle = linspace(5,100,77);

for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        error(['Parameter ' num2str(i+2) ' is not a property.']);
    end
    switch(lower(varargin{i}))
        case 'color'
            color = varargin{i+1};
        case 'threshold'
            thr = varargin{i+1};
        case 'power_norm_value'
            PowerNormValue = varargin{i+1};
        case 'smoothing'
            Smoothing = varargin{i+1};
        case 'dashed_line'
            Dashed_line = varargin{i+1};
        otherwise
            error(['Unknown property ''' num2str(varargin{i}) '''.']);
    end
end

if ~exist('color','var')
    color = 'k';
end

if ~exist('Dashed_line','var')
    Dashed_line = 1;
end

try
    PowerNormValue(PowerNormValue==0)=NaN;
end

if size(Data , 2) == 261
    Range = RangeLow;
elseif size(Data , 2) == 77
    Range = RangeMiddle;
elseif size(Data , 2) == 32
    Range = RangeHigh;
elseif size(Data , 2) ==94
    Range = RangeVHigh;
elseif size(Data , 2) ==249
    Range = RangeLow2;
else
    Range = RangeLow(end-size(Data,2)+1:end);
end

if ~exist('thr','var')
    if size(Data , 2) == 261
        thr = 13;
    elseif size(Data , 2) == 77
        thr = 21;
    elseif size(Data , 2) == 32
        thr = 13;
    elseif size(Data , 2) ==94
        thr = 50;
    elseif size(Data , 2) ==249
        thr = 1;
    else
        thr=1;
    end
end

if ~exist('PowerNormValue','var')
    Data_to_use = Data./max(Data(:,thr:end)')';
    [MaxPowerValues , Freq_Max] = max(Data(:,thr:end)');
    Freq_Max = Range(Freq_Max+thr-1);
    Freq_Max(Freq_Max==Range(thr)) = NaN;
else
    try
        Data_to_use = Data./PowerNormValue;
    catch
        Data_to_use = Data./PowerNormValue';
    end
    [MaxPowerValues , Freq_Max] = max(Data(:,thr:end)');
    Freq_Max = Range(Freq_Max+thr-1);
    Freq_Max(Freq_Max==Range(thr)) = NaN;
end

Conf_Inter = nanstd(Data_to_use)/sqrt(size(Data_to_use,1));
clear Mean_All_Sp; Mean_All_Sp=nanmean(Data_to_use);

if ~exist('Smoothing','var')
    h = shadedErrorBar(Range , Mean_All_Sp , Conf_Inter,'-k',1); hold on;
else
    h = shadedErrorBar(Range , runmean_BM(Mean_All_Sp,Smoothing) , runmean_BM(Conf_Inter,Smoothing) ,'-k',1); hold on;
end
h.mainLine.Color=color; h.patch.FaceColor=color; h.edge(1).Color=color; h.edge(2).Color=color;
xlabel('Frequency (Hz)'), ylabel('Power (a.u.)')

if Dashed_line
    [~,d]=max(nanmean(Data_to_use(:,thr:end)));
    a = vline(Range(d+thr-1),'--'); set(a,'Color',color)
end


