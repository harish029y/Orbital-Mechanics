% Define constants & orbital parameters
mu = 3.986e5; % Gravitational constant [km^3/s^2]
Re = 6378.137; % Radius of Earth [km] 
a = 7500; % semi-major axis [km]
e = 0; % eccentricity of circle
inc = 98; % inclination [deg]
omega = 0; % argument of perigee (AOP) [deg]
RAAN = 0; % Right Ascension of the Ascending Node [deg]

% defining latitude and longitude values and convert it to radians
lat = 0:10:90; % latitude values [deg]
lon = 0:10:90; % longitude values [deg]
lat_rad = lat*pi/180;
lon_rad = lon*pi/180;

% Initializing ground perturbations array as zeros
delta_h = zeros(length(lat),length(lon));

% By using for loop for defining ground perturbations for each latitude and longitude
for i = 1:length(lat)
    for j = 1:length(lon)
        r = a*(1 - e^2)/(1 + e*cos(lat_rad(i)));
        x = r*cos(lat_rad(i))*cos(lon_rad(j));
        y = r*cos(lat_rad(i))*sin(lon_rad(j));
        z = r*sin(lat_rad(i));
        
        % Rotate satellite position to Earth-fixed frame
        x_e = x*(cosd(omega)*cosd(inc) - sind(omega)*sind(inc)*cosd(lat_rad(i))) - ...
            y*(sind(omega)*cosd(inc) + cosd(omega)*sind(inc)*cosd(lat_rad(i))) + ...
            z*sind(inc)*sind(lat_rad(i));
        y_e = x*(cosd(omega)*sind(inc) + sind(omega)*cosd(inc)*cosd(lat_rad(i))) + ...
            y*(cosd(omega)*cosd(inc) - sind(omega)*sind(inc)*cosd(lat_rad(i))) - ...
            z*cosd(inc)*sind(lat_rad(i));
        z_e = x*sind(omega)*sind(lat_rad(i)) + y*cosd(omega)*sind(lat_rad(i)) + z*cosd(lat_rad(i));
        
        % Calculate ground perturbation value
        delta_h(i,j) = mu/norm([x_e,y_e,z_e]) - Re;
    end
end

% Plot ground perturbations
figure
plot(lon,delta_h)
xlabel('Longitude [deg]')
ylabel('Ground perturbation [km]')
title('Ground perturbations for variable orbital parameters')
legend(num2str(lat')) % To show especifically latitude values
