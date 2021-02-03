%% pendulum problem in Hamiltonian formulation
% H = p^2/(2*l^2) + g*l*(1-cos(theta))
% where p = l^2*dtheta/dt: 
% first order system in [theta;p]
% dtheta/dt = p/l^2
% dp/dt = -g*l*sin(theta)


%% problem parameters
M = 2; % dimension of problem 
options.model.function = 'pendulumHamiltonianROM';
l = 2;  % length (m)
g = 10; % gravity acceleration (m/s^2)
% initial condition needs to be prescribed for theta and p
u_start = [7*pi/8;0]; 

options.model.constants(1) = l;
options.model.constants(2) = g;

%% construct reduced basis
snapshots = load('pendulumHamiltonian_case1.mat');

% the snapshot matrix is constructed based on u_plot: 
% u_plot = [q(t1) q(t2) ... q(tN);
%           p(t1) p(t2) ... p(tN)];
% with size Ndim x Nt, where Ndim is the dimension of the vector of
% unknowns (= number of ODEs); which is 2 in this case
% following Pend & Mohsen, we construct the matrix
% S = Sq + i * Sp

I = sqrt(-1);
skip = 1;
S = snapshots.u_plot(1,1:skip:end) + I*snapshots.u_plot(2,1:skip:end);


% compute the economic SVD of S: S = K*Sigma*L'
% K is Ndim x Ndim, Sigma is Ndim x Ndim, L is Nt x Ndim
% in this case Ndim<<Nt
[K,Sigma,L] = svd(S);

% retain k leading singular values, k<=Ndim
k = 1;
% construct reduced basis from K, such that V'*V = I
V = K(:,1:k);
phi = real(V);
psi = imag(V);

% Jfull 
Jfull = [0 1; -1 0];
A     = [phi -psi; psi phi];
J     = A'*Jfull*A;
Aplus = J'*A'*Jfull;
options.model.A  = A;
options.model.J = J;

% get initial condition for the ROM
u_start = Aplus*u_start;
M = 2*k;

%% time integration settings
N_list  = 200; % number of time steps
t_start = 0;
t_end   = 20;

%% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2; % 1: FD, 2: AD
newton  = 2;