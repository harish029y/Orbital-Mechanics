% Gravitational constants for the Sun and the Moon
G = 6.67430e-11; 
M_sun = 1.989e30;
M_moon = 7.34767309e22;

% Initial positions of the Sun and the Moon (in meters)
r_sun = [-1.496e11, 0, 0]; 
r_moon = [3.844e8, 0, 0]; 
lat = 45;
lon = 60;

% Converting into radians
lat = lat*pi/180;
lon = lon*pi/180;

% Gravitational acceleration
a_sun = G*M_sun*r_sun/norm(r_sun)^3;
a_moon = G*M_moon*r_moon/norm(r_moon)^3;

% Gravitational potential
Usun = -G*M_sun/norm(r_sun);
Umoon = -G*M_moon/norm(r_moon);

% Total Value of gravitational potential
Utotal = Usun + Umoon;

% W.r.t Latitude and Longitude gravitational potential
[lat_grid, lon_grid] = meshgrid(-90:10:90, -180:10:180);
Utotal_grid = zeros(size(lat_grid));
for i = 1:size(lat_grid, 1)
    for j = 1:size(lat_grid, 2)
        lat = lat_grid(i, j)*pi/180;
        lon = lon_grid(i, j)*pi/180;
        r_sun = [-1.496e11*cos(lon)*cos(lat) - 1.496e11*sin(lat), -1.496e11*sin(lon), 1.496e11*cos(lat)*cos(lon)];
        r_moon = [3.844e8*cos(lon)*cos(lat) - 3.844e8*sin(lat), 3.844e8*sin(lon), 3.844e8*cos(lat)*cos(lon)];
        Usun = -G*M_sun/norm(r_sun);
        Umoon = -G*M_moon/norm(r_moon);
        Utotal_grid(i, j) = Usun + Umoon;
    end
end

figure;
contourf(lon_grid, lat_grid, Utotal_grid);
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
colorbar;
