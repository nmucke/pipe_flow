%% PDE:
% du/dt + u du/dx = nu*d2u/dx^2
% or equivalently
% du/dt + d (0.5*u^2) /dx = nu*d2u/dx^2
% discretized as
% du/dt + C*(u^2) = D * u

% see Rang, JCAM 2015, Improved traditional Rosenbrock–Wanner methods for stiff
% ODEs and DAEs


%% physical values
L     = 1; % length of domain
nu    = 0.001; % diffusion coefficient


%% discretization
M     = 400; % dimension of the model = number of interior grid points
options.model.function     = 'Burgers';
options.model.constants(1) = nu; % diffusion coefficient

x_start = 0;
x_end   = L;
x       = linspace(x_start,x_end,M+2)';  % grid points including boundary
x_in    = x(2:end-1); % position of unknowns
dx      = L/(M+1); % grid spacing, this equals diff(x)

diagonal = ones(M+1,1);
% diffusion matrix
D        = nu*(1/(dx^2))*spdiags([diagonal -2*diagonal diagonal],[0 1 2],M,M+2);
% convection matrix
% d/dx u^2 ~ ((u^2)_{i+1} - (u^2)_{i-1})/(2*dx)
% OR:
% d/dx u^2 ~ (u^2_{i+1/2} - u^2_{i-1/2})/(dx) = (0.5*(u_{i}+u_{i+1}))^2 - (0.5*(u_{i}+u_{i-1}))^2
% = D*((A*u).^2)
A = 0.5*spdiags([diagonal diagonal],[0 1],M+1,M+2);
% C = 0.5*(0.5/dx)*spdiags([-diagonal diagonal],[-1 1],M,M);
C = 0.5*(1/dx)*spdiags([-diagonal diagonal],[0 1],M,M+1);

% transfer to structure which will be available in functions
options.model.x = x_in;
options.model.x_start = x_start;
options.model.x_end = x_end;
options.model.M = M;

options.model.D = D;
options.model.A = A;
options.model.C = C;

%% initial condition
t_start = 0;
t_end   = 1;

u_start = Burgers_ex(t_start,options);

% number of time steps
N_list  = [1000];% 20 40 80];

% note: dt = (t_end - t_start)/N_list

% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2; % 1: FD, 2: AD
newton  = 2;