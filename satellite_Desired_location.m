% At the starting phrase of the satellite
altitude_initial = 850e3;
inclination_initial = input('Enter the initial input value'); % radians
% Simulation parameters0.
duration = 26060; % seconds
dt = 60; % seconds
time = 0:dt:duration;
%  Constants for the define values
mu = 3.986004418e14;
R_earth = 6.371e6;
% Calculate semi-major axis and velocity for initial altitude
a = (altitude_initial + R_earth) / (1 - cos(inclination_initial));
v_initial = sqrt(mu * (2/altitude_initial - 1/a));
position = zeros(length(time), 3);
velocity = zeros(length(time), 3);
% Fixing the initial position and velocity
position(1,:) = [altitude_initial * cos(inclination_initial), altitude_initial * sin(inclination_initial), 0];
velocity(1,:) = [-sin(inclination_initial) * v_initial, cos(inclination_initial) * v_initial, 0];
% Creating Loops through time steps
for i = 2:length(time)
% Calculate acceleration due to solar pressure
r = norm(position(i-1,:));
acceleration = (4.56e-6) / 1000 * [position(i-1,1), position(i-1,2), position(i-1,3)] / r^3;
% Updating the velocity and position from time to time
velocity(i,:) = velocity(i-1,:) + (acceleration * dt);
position(i,:) = position(i-1,:) + (velocity(i,:) * dt);
end
%Creating loops to get back the satellite from desired path