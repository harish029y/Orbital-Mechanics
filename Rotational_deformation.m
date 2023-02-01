% Initilize Parameter Values
a = 7200; % semi-major axis (km)
e = 0; % eccentricity
i = 98; % inclination (degrees)
J2 = 0.00108263;% J2 Perturbations
R = 6371; % Earth's radius (km)

% Define latitude values
lat = linspace(-90, 90, 100);

% Geopotential perturbations Equation
perturbation = (3/2)*J2*((R/a)^2)*cos(lat*pi/180).*(1 - e^2)^2*sin(i*pi/180);

% Graph Plotting
plot(lat, perturbation);
xlabel('Latitude (degrees)');
ylabel('Rotational Deformation Perturbation (m^2/s^2)');
