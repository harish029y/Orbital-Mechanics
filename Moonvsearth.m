% Define constants
G = 6.6743e-11;
M_earth = 5.9722e24;
M_moon = 7.3477e22;
M_sun = 1.989e30;

% Initial positions and velocities
r_earth = [1.4960e11; 0; 0];
v_earth = [0; 2.9783e4; 0];
r_moon = [1.4960e11+3.844e8; 0; 0];
v_moon = [0; 2.9783e4+1.022e3; 0];
r_sun = [0; 0; 0];
v_sun = [0; 0; 0];

% Defining time step and duration
dt = 3600; % time
t_max = 365*24*3600; % Maximum time in terms of duration (s)

% To store position and velocity data
t_array = 0:dt:t_max;
r_earth_array = zeros(3, length(t_array));
r_moon_array = zeros(3, length(t_array));

% Initial conditions of Earth & Moon
r_earth_array(:, 1) = r_earth;
r_moon_array(:, 1) = r_moon;

for i = 1:length(t_array)-1
    % Finding acceleration of Earth and Moon due to Sun
    a_earth_sun = -G*M_sun*(r_earth - r_sun)/norm(r_earth - r_sun)^3;
    a_moon_sun = -G*M_sun*(r_moon - r_sun)/norm(r_moon - r_sun)^3;
    
    % Finding acceleration of Earth due to Moon
    a_earth_moon = -G*M_moon*(r_earth - r_moon)/norm(r_earth - r_moon)^3;
    
    % Update velocities & Positions
    v_earth = v_earth + dt*(a_earth_sun + a_earth_moon);
    v_moon = v_moon + dt*(a_moon_sun - a_earth_moon);

    r_earth = r_earth + dt*v_earth;
    r_moon = r_moon + dt*v_moon;
    
    r_earth_array(:, i+1) = r_earth;
    r_moon_array(:, i+1) = r_moon;
end

% Finding longitude and latitude
[long, lat] = cart2sph(r_moon_array(1,:), r_moon_array(2,:), r_moon_array(3,:));

% Plot trajectory of Moon in 2D
figure
plot(r_moon_array(1,:), r_moon_array(2,:))
xlabel('x (m) -(W.r.t X-axis)')
ylabel('y (m) -(W.r.t Y-axis)')
title('Trajectory of Moon')

figure
plot(t_array/3600/24, long)
xlabel('Time (days)')
ylabel('Longitude (rad)')
title('Longitude of Moon')

% Plot latitude of Moon
figure
plot(t_array/3600/24, lat*180/pi)
xlabel('Time (days)')
ylabel('Latitude (deg)')
title('Latitude of Moon')


