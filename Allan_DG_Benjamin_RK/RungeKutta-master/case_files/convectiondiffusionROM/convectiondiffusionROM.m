%% PDE:
% du/dt = -c*du/dx + alpha*d2u/dx^2
% discretized as
% du/dt = A*u = -C*u + D*u
% with periodic boundary conditions 

%% physical values
% L     = 1; % length of domain
alpha = 0; % diffusion coefficient
c     = 1; % convection coefficient


%% discretization
M     = 25; % dimension of the model = number of interior grid points
options.model.function     = 'convectiondiffusionROM';
options.model.constants(1) = alpha; % diffusion coefficient
% options.model.constants(2) = L; % domain length
options.model.constants(3) = c; % convection coefficient

x       = linspace(-1,1,M+1)';  % grid points including boundary
x_in    = x(1:end-1);
dx      = (x(end)-x(1))/M; % grid spacing, this equals diff(x)

options.model.x_in= x_in;

% diffusion:
diagonal = ones(M,1);
D        = alpha*(1/(dx^2))*spdiags([diagonal -2*diagonal diagonal],[-1 0 1],M,M);
D(1,end) = D(1,2);
D(end,1) = D(end,end-1);

% convection:
diagonal = ones(M,1);
C        = c*(1/(2*dx))*spdiags([-diagonal diagonal],[-1 1],M,M);
C(1,end) = -C(1,2);
C(end,1) = -C(end,end-1);

A        = -C + D;
options.model.A = A;

%% initial condition
u_start = sin(pi*x_in); % initial condition

t_start = 0;
t_end   = 2.5;

N_list  = 100;

% note: dt = (t_end - t_start)/N_list

% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2; % 1: FD, 2: AD
newton  = 1; % would be nice to make a special case in the code for linear problems, where the Jacobian only has to be determined once (newton=0)


%% construct reduced basis

% skip relates to number of snapshots used; Nsnaps = Noriginal/skip
skip   = 1;
% retain k leading singular values, k<=Ndim
k      = 2;
% load snapshots from high-fidelity sim
% case1: M=25, N=100, t_end=2.5, c=1, alpha=0.1, RK4
% snapshots = load('convectiondiffusion_case1.mat');
% case2: M=25, N=100, t_end=2.5, c=0, alpha=0.1, RK4
% snapshots = load('convectiondiffusion_case2.mat');
% case3: M=25, N=100, t_end=2.5, c=1, alpha=0, RK4
snapshots = load('convectiondiffusion_case3.mat');

Morig = M; % number of nodes

S = snapshots.u_plot(:,1:skip:end);

% compute the economic SVD of S: S = K*Sigma*L'
% K is Ndim x Ndim, Sigma is Ndim x Ndim, L is Nt x Ndim
% in this case Ndim<<Nt
[K,Sigma,L] = svd(S);

% construct reduced basis from K, such that V'*V = I
V = K(:,1:k);

options.model.V = V;
options.model.Morig = Morig;

% transformed initial condition
u_start_orig = u_start;
u_start = V'*u_start;

% number of unknowns
M = k;

% precompute matrix
options.model.A = V'*A*V;

% check projection error
% proj_error = norm(speye(Morig) - V*V','fro')
proj_error2 = norm((speye(Morig) - V*V')*S,'fro')
proj_error_start = norm((speye(Morig) - V*V')*u_start_orig)

figure
semilogy(diag(Sigma)/sum(diag(Sigma)),'s');
grid
title('singular values of snapshot matrix')

% options.model.ROM_basis = ROM_basis;