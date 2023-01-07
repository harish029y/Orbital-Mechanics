% Define constants
mu1 = 3.986e5; % gravitational parameter of Earth [km^3/s^2]
mu2 = 4.90e3; % gravitational parameter of Moon [km^3/s^2]

% Initial state vector
r0 = [7000; 0]; % initial position vector [km]
v0 = [0; 8.2]; % initial velocity vector [km/s]
x0 = [r0; v0]; % initial state vector

% Time span
t0 = 0; % start time [s]
tf = 2*24*3600; % end time [s]
dt = 60; % time step [s]
t = t0:dt:tf; % time vector

% Integrate orbits using ODE45
[~,X] = ode45(@(t,x) threebodypert(t,x,mu1,mu2), t, x0);

% Extract position and velocity vectors
r = X(:,1:2); % position vectors [km]

% Plot orbit
plot(r(:,1),r(:,2));
xlabel('x'); ylabel('y');
axis equal;

% Define three-body perturbation function
function dx = threebodypert(~,x,mu1,mu2)
    % Unpack state vector
    r = x(1:2);
    v = x(3:4);
    
    % Compute accelerations
    a1 = -mu1*r/norm(r)^3; % acceleration due to Earth
    r2 = r - [1; 0]; % position vector relative to Moon
    a2 = -mu2*r2/norm(r2)^3; % acceleration due to Moon
    a = a1 + a2; % total acceleration
    
    % Pack derivative of state vector
    dx = [v; a];
end
