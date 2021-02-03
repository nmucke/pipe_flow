%% pendulum problem:
% d2 theta/ dt^2 + (g/l) sin(theta) = 0
% first order system: let u = theta, v = dtheta/dt
% du / dt = v
% dv / dt = - (g/l) sin(u)

%% problem parameters
M = 2; % dimension of problem 
options.model.function = 'pendulumROM';
l = 2;  % length (m)
g = 10; % gravity acceleration (m/s^2)
u_start = [7*pi/8;0]; % initial condition: pos (u), velocity (v)

options.model.constants(1) = l;
options.model.constants(2) = g;

%% construct reduced basis
snapshots = load('pendulum_case1.mat');

% the snapshot matrix is given by u_plot, or a subset of it
% with size Ndim x Nt, where Ndim is the dimension of the vector of
% unknowns (= number of ODEs); which is 2 in this case
skip = 1;
S = snapshots.u_plot(:,1:skip:end);

% compute the economic SVD of S: S = K*Sigma*L'
% K is Ndim x Ndim, Sigma is Ndim x Ndim, L is Nt x Ndim
% in this case Ndim<<Nt
[K,Sigma,L] = svd(S,'econ');

% retain k leading singular values, k<=Ndim
k = 2;
% construct reduced basis from K, such that V'*V = I
V = K(:,1:k);

options.model.V = V;

% get initial condition for the ROM
u_start = V'*u_start;
M = k;

%% time integration settings
N_list  = 200; % number of time steps
t_start = 0;
t_end   = 20;

%% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2; % 1: FD, 2: AD
newton  = 2;