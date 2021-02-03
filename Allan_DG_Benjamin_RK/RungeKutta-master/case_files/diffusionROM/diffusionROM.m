%% PDE:
% du/dt = nu*d2u/dx^2
% discretized as
% du/dt = D * u
% with zero Dirichlet boundary conditions at x=0 and x=L

%% physical values
L     = 10; % length of domain
alpha = 1;


%% discretization
M     = 200; % dimension of the model = number of interior grid points of full order model
options.model.function     = 'diffusionROM';
options.model.constants(1) = alpha; % diffusion coefficient
options.model.constants(2) = L; % domain length

x       = linspace(0,L,M+2)';  % grid points including boundary
x_in    = x(2:end-1);
dx      = L/(M+1); % grid spacing, this equals diff(x)

diagonal = ones(M,1);
D        = alpha*(1/(dx^2))*spdiags([diagonal -2*diagonal diagonal],[-1 0 1],M,M);

options.model.D = D;

%% initial condition
u_start = sin(pi*x_in/L); % initial condition

t_start = 0;
t_end   = 100;

N_list  = 100;

% note: dt = (t_end - t_start)/N_list

% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2; % 1: FD, 2: AD
newton  = 1;


%% construct reduced basis

% skip relates to number of snapshots used; Nsnaps = Noriginal/skip
skip   = 1;
% retain k leading singular values, k<=Ndim
k      = 1;
% load snapshots from high-fidelity sim
snapshots = load('diffusion_case1.mat');


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
options.model.D = V'*D*V;

% check projection error
% proj_error = norm(speye(Morig) - V*V','fro')
proj_error2 = norm((speye(Morig) - V*V')*S,'fro')
proj_error_start = norm((speye(Morig) - V*V')*u_start_orig)


figure
semilogy(diag(Sigma)/sum(diag(Sigma)),'s');
grid
title('singular values of snapshot matrix')

% options.model.ROM_basis = ROM_basis;