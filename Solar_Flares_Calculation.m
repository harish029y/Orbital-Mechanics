clc
clear

% Defining initial Values
altitude_initial = 850e3; % m
inclination_initial = input('Enter the initial input value'); % rad
mu = 3.986004418e14; % m^3/s^2
R_earth = 6.371e6; % m

% Time initilization parameters
duration = 26060; % s
dt = 60; % s
time = 0:dt:duration;

% Calculating semi-major axis from definite altitudes
a = (altitude_initial + R_earth) / (1 - cos(inclination_initial));

% Arrays to store position and velocity for the definite period
position = zeros(length(time), 3);
velocity = zeros(length(time), 3);

position(1,:) = [altitude_initial * cos(inclination_initial), altitude_initial * sin(inclination_initial), 0];

% Finding the initial velocity at the initial altitude
v_initial = sqrt(mu * (2/norm(position(1,:)) - 1/a));

% Set initial position and velocity
position(1,:) = [altitude_initial * cos(inclination_initial), altitude_initial * sin(inclination_initial), 0];
velocity(1,:) = [-sin(inclination_initial) * v_initial, cos(inclination_initial) * v_initial, 0];

% For Solar Pressure perturbations, so we use gravitational attraction values
mu_sun = 1.32712440018e20; % m^3/s^2
R_sun = 6.96e8; % m

% Initialize position vector's to store array values
position_vector = zeros(length(time), 3);
position_vector(1,:) = position(1,:);

% Solar Flare Parameters
solar_flare_probability = 0.05; % Probability of solar flare occurrence
max_solar_flare_magnitude = 0.15; % Maximum solar flare magnitude

% For loop to calculate the effects of solar flares on the satellite's attitude
for i = 2:length(time)
    % Acceleration due to solar pressure
    r = norm(position(i-1,:));
    acceleration_sp = (4.56e-6) / 1000 * [position(i-1,1), position(i-1,2), position(i-1,3)] / r^3;
    
    % Acceleration due to sun's gravitational attraction
    r_vec_to_sun = -position(i-1,:);
    r_to_sun = norm(r_vec_to_sun);
    acceleration_grav = -mu_sun * r_vec_to_sun / r_to_sun^3;
    
    % Total acceleration will be always calculated with the help of addition from Solar Pressure & gravitational attraction
    acceleration = acceleration_sp + acceleration_grav;
    
% Check for the occurrence of a solar flare
    if rand(1) <= solar_flare_probability
        solar_flare_magnitude = max_solar_flare_magnitude * rand(1);
        acceleration = acceleration + solar_flare_magnitude * acceleration;
        
        % Modify the satellite attitude due to a solar flare
        perturbation_vector = rand(1,3) - 0.5; % Random perturbation vector
        attitude_perturbation = cross(velocity(i,:), perturbation_vector) / norm(velocity(i,:)) * solar_flare_magnitude * dt;
        position(i,:) = position(i,:) + attitude_perturbation;
        % Recalculate the inclination angle
        inclination_angle = atan2(position(i,2), position(i,1));
    end
end
for i = 2:length(time)
    % Acceleration due to solar pressure
    r = norm(position(i-1,:));
    acceleration_sp = (4.56e-6) / 1000 * [position(i-1,1), position(i-1,2), position(i-1,3)] / r^3;
    
    % Acceleration due to sun's gravitational attraction
    r_vec_to_sun = -position(i-1,:);
    r_to_sun = norm(r_vec_to_sun);
    acceleration_grav = -mu_sun * r_vec_to_sun / r_to_sun^3;
    
    % Total acceleration will be always calculated with the help of addition from Solar Pressure & gravitational attraction
    acceleration = acceleration_sp + acceleration_grav;
    
    % Check for the occurrence of a solar flare
    if rand(1) <= solar_flare_probability
        solar_flare_magnitude = max_solar_flare_magnitude * rand(1);
        acceleration = acceleration + solar_flare_magnitude * acceleration;
%{
        % Modify the satellite attitude due to a solar flare
        perturbation_vector = rand(1,3) - 0.5; % Random perturbation vector
        attitude_perturbation = cross(position(i,:), perturbation_vector) / norm(position(i,:)) * solar_flare_magnitude * dt;
        position(i,:) = position(i,:) + attitude_perturbation;
%}
        % Modify the satellite position due to a solar flare
         perturbation_vector = rand(1,3) - 0.5; % Random perturbation vector
         position(i,:) = position(i,:) + perturbation_vector * solar_flare_magnitude * dt;
    
    end
    
    % Update position and velocity
    position(i,:) = position(i-1,:) + velocity(i-1,:) * dt + 0.5 * acceleration * dt^2;
    velocity(i,:) = velocity(i-1,:) + acceleration * dt;
end


% Plotting the relationship between inclination and solar flares
inclination = zeros(length(time), 1);
solar_flare = zeros(length(time), 1);
for i = 1:length(time)
    inclination(i) = atan2(position(i,2), position(i,1));
    if rand(1) <= solar_flare_probability
        solar_flare(i) = max_solar_flare_magnitude * rand(1);
    end
end


% Plotting the relationship between inclination and solar flares
inclination = zeros(length(time), 1);
solar_flare = zeros(length(time), 1);
for i = 1:length(time)
inclination(i) = atan2(position(i,2), position(i,1));
if rand(1) <= solar_flare_probability
solar_flare(i) = 1;
end
end

figure()
subplot(2,1,1)
plot(time/3600, inclination*180/pi)
xlabel('Time (hours)')
ylabel('Inclination (deg)')
title('Relationship between Inclination and Solar Flares')
subplot(2,1,2)
plot(time/3600, solar_flare)
xlabel('Time (hours)')
ylabel('Solar Flare')
title('Occurrence of Solar Flares over Time')

% Plotting the relationship between inclination and solar flares
figure
plot(inclination, solar_flare, '.', 'MarkerSize', 10);
xlabel('Inclination (rad)');
ylabel('Solar Flare Magnitude (m/s^2)');
title('Relationship between Inclination and Solar Flares');

% Plotting the relationship between inclination and solar flares
inclination = zeros(length(time), 1);
solar_flare = zeros(length(time), 1);
for i = 1:length(time)
inclination(i) = atan2(position(i,2), position(i,1));
if rand(1) <= solar_flare_probability
solar_flare(i) = max_solar_flare_magnitude * rand(1);
end
end

% Create scatter plot
scatter(inclination, solar_flare, 'filled')
xlabel('Inclination')
ylabel('Solar Flare Magnitude')

%Plotting the orbit of the satellite
figure
% First obtain the radius of the Earth
earth_radius = 6.371e6; % m
% Set the aspect ratio for the plot
% The x, y, and z axis will have the same scale.
daspect([1 1 1])

% Define a sphere to represent the Earth
[X, Y, Z] = sphere(50);

% Scaling down the size of the earth to provide better presentation
X = earth_radius / 10^3 * X;
Y = earth_radius / 10^3 * Y;
Z = earth_radius / 10^3 * Z;

% Display the Earth
surf(X, Y, Z, 'EdgeColor', 'none')
hold on

% Plot the orbit of the satellite
plot3(position(:,1), position(:,2), position(:,3), 'r', 'LineWidth', 2)

% Add title and labels
title('Orbit of the satellite')
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')

% Set axis limits to include the whole orbit
axis([-2 2 -2 2 -2 2] * (a + earth_radius) / 10^3)

% Add a legend
legend('Orbit')