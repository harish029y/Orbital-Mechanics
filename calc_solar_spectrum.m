% Defining constants related to the spacecraft
AU = 1.495978707e11; % Iniate AU
r_planet = 6.371e6; % Earth Radius
r_sun = 6.957e8; % Sun Radius
d = 1.5*AU; % Earth to Sun distance in AU
lambda_min = 400e-9;
lambda_max = 700e-9;
n_lambda = 100;
lambda = linspace(lambda_min, lambda_max, n_lambda);

% Solar Spectrum Calculation
I_0 = solar_spectrum(lambda);

% Earth Perturbation Calculation
I = zeros(size(lambda));
for i = 1:n_lambda
  I(i) = I_0(i)*exp(-2*r_planet/AU*(r_sun/d)^2/lambda(i));
end

% plot Graph
figure
plot(lambda, I, 'r', 'LineWidth', 2)
xlabel('Wavelength (m)')
ylabel('Intensity (W/m^2/sr/m)')

function I_0 = solar_spectrum(lambda)
  h = 6.62607004e-34; % Planck constant
  c = 299792458; % speed of light
  k = 1.38064852e-23; % Boltzmann constant
  T = 5778; % Solar Temperature
  
  % Solar spectrum Calculation
  I_0 = 2*h*c^2./lambda.^5./(exp(h*c./(lambda*k*T)) - 1);
end
