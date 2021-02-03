%% pendulum problem:
% d2 theta/ dt^2 + (g/l) sin(theta) = 0
% first order system: let u = theta, v = dtheta/dt
% du / dt = v
% dv / dt = - (g/l) sin(u)

%% problem parameters
M = 2; % dimension of problem 
options.model.function = 'pendulum';
l = 2;  % length (m)
g = 10; % gravity acceleration (m/s^2)
u_start = [7*pi/8;0]; % initial condition: pos (u), velocity (v)

options.model.constants(1) = l;
options.model.constants(2) = g;


%% time integration settings
N_list  = 200; % number of time steps
t_start = 0;
t_end   = 20;

%% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2; % 1: FD, 2: AD
newton  = 2;