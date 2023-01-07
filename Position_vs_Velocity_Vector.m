% Create random values for latitude & longitude data, random position & velocity vectors
lat = rand(100, 1) * 180 - 90; % Latitude values ranging from -90 to 90
lon = rand(100, 1) * 360 - 180; % Longitude values ranging from -180 to 180
pos = rand(100, 3); % position vectors
vel = rand(100, 3); % velocity vectors

% Specify the orbital parameters
a = 7500; % Semimajor axis
e = 0; % Eccentricity of circle
i = 98; % Inclination of angle [deg]
RAAN = 0; % Right ascension of the ascending node [deg]
w = 0; % Argument of perigee [deg]
nu = 90; % True anomaly
AOP = 0; % Argument of periapsis [deg]

% 2D plot - latitude and longitude data
figure;
plot(lon, lat, '.')
xlabel('Longitude')
ylabel('Latitude')
title('Latitude vs Longitude')

% Updating position and velocity vectors to the plot
hold on
quiver(lon, lat, pos(:,1), pos(:,2))
quiver(lon, lat, vel(:,1), vel(:,2), 'r')
hold off
legend('Data points', 'Position vectors', 'Velocity vectors')

% Displaying Orbital parameters along with their values
disp('Orbital parameters:')
disp(['Semimajor axis: ', num2str(a), ' km'])
disp(['Eccentricity: ', num2str(e)])
disp(['Inclination: ', num2str(i), ' degrees'])
disp(['Right ascension of the ascending node: ', num2str(RAAN), ' degrees'])
disp(['Argument of perigee: ', num2str(w), ' degrees'])
disp(['True anomaly: ', num2str(nu), ' degrees'])

% Uploading latitude and longitude data
lat = rand(100, 1) * 180 - 90; % Latitude values ranging from -90 to 90
lon = rand(100, 1) * 360 - 180; % Longitude values ranging from -180 to 180

% Plotting 2D plot of the latitude and longitude data
figure;
plot(lon, lat, '.')
xlabel('Longitude')
ylabel('Latitude')
title('Latitude vs Longitude')
