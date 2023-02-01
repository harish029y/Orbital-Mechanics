% Initilize Parameter Values
a = 7200; % semi-major axis (km)
e = 0; % eccentricity
i = 98; % inclination (degrees)
k2 = 0.299;
R = 6371; % Radius of Earth
G = 6.67*10^-11; % gravitational constant
M_E = 5.97*10^24; % Earth's mass (kg)
M_m = 7.35*10^22; % Moon's mass (kg)
D = 384400; % distance between Moon & Earth (km)

% Define latitude values
lat = linspace(-90, 90, 100);

% Ocean tide perturbations Equation
perturbation = (3/2)*k2*((R/a)^3)*(G*M_E/D^3)*(1-3*sin(lat*pi/180).^2)*sin(i*pi/180);

% Plot the perturbations
plot(lat, perturbation);
xlabel('Latitude (degrees)');
ylabel('Ocean Tide Perturbation (m^2/s^2)');
