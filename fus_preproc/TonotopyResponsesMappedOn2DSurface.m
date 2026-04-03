%% Step 1: load the 4 tonotopy slices with the wanted statistics to project
% (ex; best frequency response)
datapath = sessions{1};
MatBestFreq = load(fullfile(datapath, 'MatBestFreq.mat'));
MatBestFreq = MatBestFreq.MatBestFreq;
%% For each slice
Slices = {'A', 'B', 'C', 'D'};
Results = struct();
for iSlice = 1:length(Slices)
% 1: Draw cortical ROI 
    fUSImage = MatBestFreq.(Slices{iSlice}).MeanImage;
    
    % Log scale 
    frame_double = double(fUSImage);
    frame_log = log(frame_double + 1);
    frame_log = frame_log - min(frame_log(:));
    frame_log = frame_log / max(frame_log(:));
    frame_log = uint8(frame_log * 255);
    Results.(Slices{iSlice}).MeanImageLog = frame_log;
    figure(1);
    imagesc(frame_log)
    axis image
    colormap hot
    title('Draw broad cortical mask');
    h = drawpolygon;
    Results.(Slices{iSlice}).ACmask = createMask(h);
end


%% 2: extract cortical top layer with position of top voxels
figure('Units', 'normalized', 'Position', [0.1, 0.5, 0.7, 0.3]);
for iSlice = 1:length(Slices)
    mask = Results.(Slices{iSlice}).ACmask;
    [y, x] = find(mask);

    % For each column, keep highest (smallest y)
    unique_x = unique(x);
    top_points = [];

    for i = 1:length(unique_x)
        col = unique_x(i);
        ys = y(x == col);
        top_y = min(ys); % top of cortex
        top_points = [top_points; col, top_y];
    end
    
    Results.(Slices{iSlice}).top_points = top_points;
    X = top_points(:,1);
    Y = top_points(:,2);
    p = polyfit(X, Y, 1); % linear fit
    slope = p(1);
    theta = atan(slope); % rotation angle (radians)
    Results.(Slices{iSlice}).theta = theta;
    subplot(1,4,iSlice)
    img = MatBestFreq.(Slices{iSlice}).MeanImage;
    imagesc(img); 
    hold on; 
    x_fit = 1:size(img,2);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'b-', 'LineWidth', 2);
    plot(X, Y, 'g.', 'MarkerSize', 2);
    
end

%% 4: rotate image 

figure('Units', 'normalized', 'Position', [0.1, 0.5, 0.7, 0.3]);

for iSlice = 1:length(Slices)
    mask = Results.(Slices{iSlice}).ACmask;
    theta = Results.(Slices{iSlice}).theta;
    img = MatBestFreq.(Slices{iSlice}).MeanImage;
    top_points = Results.(Slices{iSlice}).top_points;
    X = top_points(:,1);
    Y = top_points(:,2);
    
    
    theta_deg = rad2deg(theta);
    rot_img = imrotate(img, theta_deg, 'bilinear', 'crop');
    rot_Tonotopy = imrotate(MatBestFreq.(Slices{iSlice}).BestFrequencyMap, theta_deg, 'bilinear', 'crop');
    rot_mask = imrotate(mask, theta_deg, 'nearest', 'crop');

    % Center of image
    cx = (size(img,2)+1)/2;
    cy = (size(img,1)+1)/2;

    % Center coordinates
    Xc = X - cx;
    Yc = Y - cy;

    % Rotation matrix (NOW +theta)
    R = [cos(-theta) -sin(-theta); 
         sin(-theta)  cos(-theta)];

    rot_coords = R * [Xc'; Yc'];

    % Shift back
    Xr = rot_coords(1,:)' + cx;
    Yr = rot_coords(2,:)' + cy;
    
    
    % Fit check
    p_rot = polyfit(Xr, Yr, 1);
    x_fit = linspace(min(Xr), max(Xr), 200);
    y_fit = polyval(p_rot, x_fit);

    
    ax1= subplot(4,4,iSlice);
    imagesc(ax1, img);
    axis image;
    colormap(ax1, 'hot');
        
    ax3 = subplot(4,4,iSlice+4);
    imagesc(rot_img);
    axis image;
    colormap(ax3, 'hot');
    hold on;

    [y_m, x_m] = find(rot_mask);
    plot(x_m, y_m, 'g.', 'MarkerSize', 5);

    plot(Xr, Yr, 'b.', 'MarkerSize', 10);
    
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
    
    
    Results.(Slices{iSlice}).RotatedMask = rot_mask;
    Results.(Slices{iSlice}).RotatedImage = rot_Tonotopy;
    Results.(Slices{iSlice}).RotatedTonotopy = rot_img;
    Results.(Slices{iSlice}).p_rot = p_rot;
    
    
    % Projection
    ImageToProject = rot_Tonotopy;
    [y_m, x_m] = find(rot_mask);
    vals = ImageToProject(sub2ind(size(ImageToProject), y_m, x_m));
    slope = p_rot(1);
    theta_rot = atan(slope);

    u_t = [cos(theta_rot); sin(theta_rot)];   % along line
    
    % Center (important!)
    cx_r = (size(ImageToProject,2)+1)/2;
    cy_r = (size(ImageToProject,1)+1)/2;

    coords = [x_m - cx_r, y_m - cy_r];

    t = coords * u_t;   % position along regression line
    nbins = length(X);
    t_edges = linspace(min(t), max(t), nbins+1);

    proj_map = nan(1, nbins);
    for b = 1:nbins
        idx = t >= t_edges(b) & t < t_edges(b+1);
        
        if any(idx)
            proj_map(b) = mean(vals(idx));
        end
    end
    
    ax2= subplot(4,4,iSlice+8);
    imagesc(ax2, ImageToProject)
    FREQSAll =[200 400 800 1600 3200 6400 12800 25600];

    cmapTonotopy = hsv(9);
    cmapTonotopy(1,:) = [0 0 0];
    nFAll = numel(FREQSAll);
    colormap(ax2, cmapTonotopy);
    caxis(ax2, [0 nFAll]);
    colorbar(ax2, 'Ticks', 0:nFAll, 'TickLabels', string([0, FREQSAll]))
    
    
    
    subplot(4,4,iSlice+12)
    plot(proj_map, 'LineWidth', 2);
    title(['Projection - ' Slices{iSlice}]);
    xlabel('Cortical position');
    ylabel('Best Frequency');
    
     Results.(Slices{iSlice}).Projection = [X, proj_map'];
end


%% Vizualisation


Nx = size(MatBestFreq.A.BestFrequencyMap, 2);
FinalMap = NaN(Nx, length(Slices));
for iSlice=1:length(Slices)
    ProjCoordinates = Results.(Slices{iSlice}).Projection;
    FinalMap(ProjCoordinates(:, 1), iSlice) = ProjCoordinates(:,2);
    
end

figure;
imagesc(FinalMap); %1:Nx, 1:size(FinalMap,1),
axis xy;          % IMPORTANT → puts origin at bottom-left
ylabel('x( pixel)');
xlabel('Slice');
colormap(cmapTonotopy);
caxis( [0 nFAll]);
colorbar('Ticks', 0:nFAll, 'TickLabels', string([0, FREQSAll]))

    
xticks(1:4)
xticklabels(Slices);
colorbar;

saveas(gcf, fullfile(datapath, 'figures/CorticalTonotopyMappedOn2DSurface.png'))

 
