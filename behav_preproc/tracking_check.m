function tracking_check(datapath, marker, range)
% range in seconds
% marker as in DLC_data.mat like  'pupil_center_007_mvt'

load(fullfile(datapath,'video','DLC_data.mat'),marker);
[~, sessionName ] = fileparts(datapath);
%  pre-existing data
marker_time = Range(eval(marker),'s'); 
marker_data = Data(eval(marker));     
marker_data = marker_data(:, 1);
videoT0 = marker_time(1);                     % global time of first video frame

if contains(marker, '007')
    camFile = dir(fullfile(datapath,'video','*007*DLC*.mp4'));
else
    camFile = dir(fullfile(datapath,'video','*004*DLC*.mp4'));
end
vr = VideoReader(fullfile(camFile.folder,camFile.name));
outFile = fullfile(datapath,'video',[marker '_tracking_check.mp4']);

range_sample = range .* vr.FrameRate; 

vw = VideoWriter(outFile,'MPEG-4');  vw.FrameRate = vr.FrameRate;  open(vw);

fig = figure('Color','w','Units','pixels','Position',[1 1 1920 1080], 'visible', 'on');
sgtitle([marker ', ' sessionName]) 
axTar = subplot(3,4,[2 3 ; 6 7]);
axR   = subplot(3,4,9:12);

% iterate over pupil samples and show the matching frame
win   = range_sample(1):range_sample(2);               % ±100 samples for the trace
% win   = -100:100;               % ±100 samples for the trace

for k = range_sample(1):range_sample(2)
    tg = marker_time(k);           % current sample in global time
    tv = tg - videoT0;      % same moment in the video clock
    
    % grab and display the corresponding frame
    vr.CurrentTime = tv;
    imshow(readFrame(vr),'Parent',axTar);  axis(axTar,'off');
    
    % draw 200-sample pupil trace centred on 'k'
    plot(axR, marker_time(win), marker_data(win));  hold(axR,'on');
    %     idx = k+win;
%     plot(axR, marker_time(idx), marker_data(idx));  hold(axR,'on');
    xline(axR, tg, 'LineWidth',2);    % vertical bar at the current time
    hold(axR,'off');
    
    % optionally capture the composite frame for a new movie
    frame = getframe(fig);
%     writeVideo(vw, frame);
end

close(vw);
end