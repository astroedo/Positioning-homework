%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Positioning and Location Based Services
% A.A. 2024/2025
% Exercise 3:  Ionospheric delay
% 
% Edoardo Pessina
% 
% References:
%     - http://www.navipedia.net/index.php/Clock_Modelling
%     - http://www.ajgeomatics.com/SatelliteOrbitsRev0828-2010.pdf
%     - http://www.navipedia.net/index.php/GPS_and_Galileo_Satellite_Coordinates_Computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
%% Import data 

%load parameters
line1 = 'GPSA 7.4506D-09 1.4901D-08 -5.9605D-08 -1.1921D-07 IONOSPHERIC CORR';
line2 = 'GPSB 9.2160D+04 1.3107D+05 -6.5536D+04 -5.2429D+05 IONOSPHERIC CORR';
ionoparams = [cell2mat(textscan(line1, '%*s %f %f %f %f %*s')) ...
cell2mat(textscan(line2, '%*s %f %f %f %f %*s'))];

%% Task 3.1 Zenithal maps of Iono corrections

%initialize values for the zenital cycle
t = [0,6*3600,12*3600,18*3600]; % [seconds]

elev = 90; % [degree]
long = [(-180):0.5:(180)]; % [degree]
lat = [-80:0.5:80]; % [degree]

%initialize matrix
delay = zeros(length(lat),length(long));

%time, phi and lambda cycle
for i = 1:length(t)
    for j = 1: length(lat)
        for k = 1:length(long)
            delay(j,k,i) = iono_correction(lat(j), long(k), 0, elev, t(i), ionoparams);
        end
    end
end

%% plots
% Create a figure with 2x2 subplots
figure('Position', [100, 100, 1600, 800]);

minDelay = min(delay(:));
maxDelay = max(delay(:));

load coastlines

for i = 1:4
    subplot(2,2,i)
    
    % First plot the coastlines
    plot(coastlon, coastlat, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    hold on;
    
    % Create meshgrid
    [LON, LAT] = meshgrid(long, lat);
    
    % Plot delay data with pcolor and some transparency
    h = pcolor(LON, LAT, delay(:,:,i));
    set(h, 'FaceAlpha', 0.7);  % Make the delay plot somewhat transparent
    shading interp;
    
    % Set consistent color limits
    caxis([minDelay maxDelay]);
    
    % Add colorbar
    colorbar;
    colormap(jet);
    
    % Set proper axis limits and labels
    axis([-180 180 -80 80]);
    xlabel('Longitude (degrees)');
    ylabel('Latitude (degrees)');
    title(sprintf('Time = %d hours', t(i)/3600));
    
    % Add grid
    grid on;
    
    % Keep the coastlines visible
    set(gca, 'Layer', 'top');
    
    hold off;
end

sgtitle('Ionosphere delay during the day - Latitude vs Longitude');

%% Task 3.2

% Milano position in degrees
phi2 = 45 + 28 / 60 + 38.28 / 60^2; %degrees
lambda2 = 9 + 10 / 60 + 53.40 / 60^2; %degrees


%initialize values for the zenital cycle
t = [0,12*3600]; % [seconds]

elev = [0:0.5:90];% [degree]
AZ = [-180:0.5:180]; % [degree]

%initialize matrix
delay2 = zeros(length(elev),length(AZ));

%time, phi and lambda cycle
for i = 1:length(t)
    for j = 1: length(AZ)
        for k = 1:length(elev)
            delay2(k,j,i) = iono_correction(lambda2, phi2, AZ(j), elev(k), t(i), ionoparams);
        end
    end
end

%% plot 2

% Create figure
figure('Position', [100, 100, 1200, 500]);

% Convert to radians for calculation
az_rad = deg2rad(AZ);
[AZ_grid, ELEV_grid] = meshgrid(AZ, elev);

% Calculate x and y coordinates for cartesian plotting
% 90-ELEV_grid to have 90° in center and 0° at edges
X = (90-ELEV_grid) .* cosd(AZ_grid);
Y = (90-ELEV_grid) .* sind(AZ_grid);

% Plot for each time (0 hours and 12 hours)
for i = 1:2
    subplot(1,2,i)
    
    % Find min and max for each time separately
    minDelay = min(min(delay2(:,:,i)));
    maxDelay = max(max(delay2(:,:,i)));
    
    % Plot using pcolor
    h = pcolor(X, Y, delay2(:,:,i));
    shading interp;
    
    alpha(h, 0.7); 

    % Make plot circular
    axis equal;
    axis([-90 90 -90 90]);
    
    % Add elevation circles
    hold on;
    theta = linspace(0, 2*pi, 100);
    for r = [30 60]
        plot(r*cos(theta), r*sin(theta), 'k:', 'LineWidth', 0.5);
        % Add elevation labels
        text(0, r, [num2str(90-r) '°'], 'HorizontalAlignment', 'left');
    end
    
    % Add azimuth labels
    r = 95;
    for az = [0 90 180 270]
        text(r*cosd(az), r*sind(az), [num2str(az) '°'], 'HorizontalAlignment', 'center');
    end
    
    % Set colorbar and limits for each plot separately
    c = colorbar;
    ylabel(c, 'Delay [meters]');
    colormap(jet);  

    caxis([minDelay maxDelay]);
    
    % Add title with time and delay range
    title(sprintf('Milano Ionospheric Delay at t = %d hours\n', ...
        t(i)/3600));
    
    % Add grid
    grid on;
    
    hold off;
end

% Add overall title
sgtitle('Milano Ionospheric Delay - Elevation vs Azimuth');