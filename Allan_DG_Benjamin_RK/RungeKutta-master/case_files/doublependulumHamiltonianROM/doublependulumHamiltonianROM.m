%% double pendulum problem in Hamiltonian formulation

% see e.g. https://diego.assencio.com/?index=e5ac36fcb129ce95a61f8e8ce0572dbf
% or https://mse.redwoods.edu/darnold/math55/DEProj/Sp00/FranScott/finalpaper.pdf
% or http://scienceworld.wolfram.com/physics/DoublePendulum.html
% or http://iopscience.iop.org/article/10.1088/1742-6596/739/1/012066/pdf
% dtheta1/dt = (L2*p1 - L1*p2*cos(theta1-theta2))/(L1*h3)
% dtheta2/dt =( -m2*L2*p1*cos(theta1-theta2) + (m1+m2)*L1*p2 )/(m2*L2*h3)
% dp1/dt = -(m1+m2)*g*l1*sin(theta1) - h1 + h2
% dp2/dt = -m2*g*l2*sin(theta2) + h1 - h2
% h1 = (p1*p2*sin(theta1-theta2))/h3
% h2 =
% sin(2*(theta1-theta2))*(m2*(L2^2)*(p1^2)+(m1+m2)*(L1^2)*(p2^2)-2*m2*L1*L2*p1*p2*cos(theta1-theta2))/
% (2*L1*L2*h3) % there is a factor 2  missing in the second and third reference?
% h3 = L1*L2*(m1+m2*sin(theta1-theta2)^2)

% p1 = (m1+m2)*(L1^2)*(dtheta1/dt) + m2*L1*L2*(dtheta2/dt)*cos(theta1-theta2)
% p2 = m2*(L2^2)*(dtheta2/dt) + m2*L1*L2*(dtheta1/dt)*cos(theta1-theta2)

%% problem parameters
M = 4; % dimension of problem 
options.model.function = 'doublependulumHamiltonianROM';
L1 = 2;  % length (m)
L2 = 1;
m1 = 1;
m2 = 1;
g = 10; % gravity acceleration (m/s^2)
% initial condition for theta and p
% set initial velocity to zero, so p1 and p2 are zero
u_start = [pi/2;-pi/2;0;0]; 

options.model.constants(1) = L1;
options.model.constants(2) = L2;
options.model.constants(3) = m1;
options.model.constants(4) = m2;
options.model.constants(5) = g;

%% construct reduced basis
snapshots = load('doublependulumHamiltonian_case1.mat');

% the snapshot matrix is constructed based on u_plot: 
% u_plot = [q(t1) q(t2) ... q(tN);
%           p(t1) p(t2) ... p(tN)];
% with size Ndim x Nt, where Ndim is the dimension of the vector of
% unknowns (= number of ODEs); which is 2 in this case
% following Pend & Mohsen, we construct the matrix
% S = Sq + i * Sp

I = sqrt(-1);
skip = 1;
indx_q = 1:M/2;
indx_p = M/2+1:M;
S = snapshots.u_plot(indx_q,1:skip:end) + I*snapshots.u_plot(indx_p,1:skip:end);


% compute the economic SVD of S: S = K*Sigma*L'
% K is Ndim x Ndim, Sigma is Ndim x Ndim, L is Nt x Ndim
% in this case Ndim<<Nt
[K,Sigma,L] = svd(S);

% retain k leading singular values, k<=Ndim
k = 2;
% construct reduced basis from K, such that V'*V = I
V = K(:,1:k);
phi = real(V);
psi = imag(V);

% Jfull 
Z  = zeros(2,2);
E  = eye(2);
Jfull = [Z E; -E Z]; % skew-symmetric matrix in Hamiltonian formulation
A     = [phi -psi; psi phi];
J     = A'*Jfull*A;  % 
Aplus = J'*A'*Jfull; % the symplectic inverse
options.model.A  = A;
options.model.J = J;

% get initial condition for the ROM
u_start = Aplus*u_start;
M = 2*k;

%% time integration settings
N_list  = 100; % number of time steps
t_start = 0;
t_end   = 0.5;

%% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2; % 1: FD, 2: AD
newton  = 2;