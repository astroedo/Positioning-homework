%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Positioning and Location Based Services
% A.A. 2024/2025
% Exercise 2:  GPS orbits computation
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
%% Import data from Almanac of satellite SVN 63, PRN 01 (Block IIR) for 2016, 11, 28

dt0 = -7.661711424589D-05;
dt1 = -3.183231456205D-12;
dt2 =  0.000000000000D+00;      % Clock offset
a = 5.153650835037D+03^2;       % [meters]
e = 3.841053112410D-03;         % [-]
b = a * sqrt(1 - e^2);          % [meters]
M0 = 1.295004883409D+00;        % [radians]
Omega0 = -2.241692424630D-01;   % [radians]
Omegadot = -8.386063598924D-09; % [radians/sec]
i0 = 9.634782624741D-01;        % [radians]
idot = -7.286017777600D-11;     % [radians/sec]
w0 = 9.419793734505D-01;        % [radians]
wdot = 0.0;                     % [radians/sec]

GMe = 3.986005D+14;             % [m3/s2]
OmegaEdot = 7.2921151467D-05;   % [radians]

%% Compute positions in ITRF, X, Y, Z
t_0 = 0;                         % [second] 
t_step = 60;                     % [second]
t_end = 86400;                   % [second]
n = sqrt(GMe / a^3);             % [1/s]

% Initializing vector of clock offsets and positions.
num_steps = floor(t_end / t_step) + 1;
X = zeros(1, num_steps);
Y = zeros(1, num_steps);
Z = zeros(1, num_steps);
Offset = zeros(1, num_steps);

for idx = 1:num_steps
    Dt = (idx - 1) * t_step;
    
    % Compute the clock offset
    Offset(idx) = dt0 + dt1 * Dt + dt2 * Dt^2;
    
    % Compute M(t)
    M_t = M0 + n * Dt;
    
    % Solve Kepler's equation iteratively
    Eta_t = M_t;
    for i = 1:10
        Eta_t = M_t + e * sin(Eta_t);
    end
    
    % Compute true anomaly
    psi_t = atan2(sqrt(1 - e^2) * sin(Eta_t), cos(Eta_t) - e);
    
    % Compute r
    r_t = a * (1 - e^2) / (1 + e * cos(psi_t));
    
    % Compute position in orbital plane
    x_Dt = r_t * cos(psi_t);
    y_Dt = r_t * sin(psi_t);
    
    % Compute w(t), i(t), and Ω(t)
    w_t = w0 + wdot * Dt;
    i_t = i0 + idot * Dt;
    W_t = Omega0 + (Omegadot - OmegaEdot) * Dt - OmegaEdot * t_0;
    
    % Define the three rotation matrices
    [R1, R3_w, R3_W] = ORStoITRF(-i_t, -w_t, -W_t);
    Coord = R3_W * R1 * R3_w * [x_Dt; y_Dt; 0];
    
    X(idx) = Coord(1);
    Y(idx) = Coord(2);
    Z(idx) = Coord(3);
end

%% Transform the set of Cartesian coordinates to a set of Geodetic coordinates
Position = [X; Y; Z];
Pos_geod = CarttoGeod(Position, a, e, b);
longitude = rad2deg(Pos_geod(1,:));
latitude = rad2deg(Pos_geod(2,:));
height = Pos_geod(3,:) ./ 1000;
height_mean = mean(height);

%% Plot geodetic coordinates on a world map
figure(1);
worldmap('World');
load('coastlines', 'coastlat', 'coastlon');
plotm(coastlat, coastlon, 'k');
plotm(latitude, longitude, 'r', 'LineWidth', 2);
title('Satellite Ground Track %f', height_mean);
xlabel('Longitude');
ylabel('Latitude');

%% Plot the height oscillation
figure(2);
plot((0:t_step:t_end)/3600, height);
xlabel('Time (hours)', 'FontWeight', 'bold');
ylabel('Height (km)', 'FontWeight', 'bold');
title(sprintf('Satellite Height Oscillation Over Time (Mean Height: %.2f km)', height_mean), 'FontWeight', 'bold');
grid on;

%% Plot the satellite clock offset
figure(3);
plot((0:t_step:t_end)/60, Offset);
xlabel('Time (minutes)', 'FontWeight', 'bold');
ylabel('Clock Offset (seconds)', 'FontWeight', 'bold');
title('Satellite Clock Offset Over Time', 'FontWeight', 'bold');
grid on;
