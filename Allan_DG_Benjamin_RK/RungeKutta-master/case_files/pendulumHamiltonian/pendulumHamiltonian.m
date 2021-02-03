%% pendulum problem in Hamiltonian formulation
% H = p^2/(2*l^2) + g*l*(1-cos(theta))
% where p = l^2*dtheta/dt: 
% first order system in [theta;p]
% dtheta/dt = p/l^2
% dp/dt = -g*l*sin(theta)


%% problem parameters
M = 2; % dimension of problem 
options.model.function = 'pendulumHamiltonian';
l = 2;  % length (m)
g = 10; % gravity acceleration (m/s^2)
% initial condition needs to be prescribed for theta and p
u_start = [7*pi/8;0]; 

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