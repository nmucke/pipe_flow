%% PDE:
% du/dt = alpha*d2u/dx^2 + f(x,t)
% discretized as
% du/dt = D * u + F
% with time-dependent Dirichlet boundary conditions at x=0 and x=L
% see Borggaard & Iliescu, ﻿Approximate deconvolution boundary conditions
% for large eddy simulation, 2006

%% physical values
L     = 1; % length of domain
alpha = 1;


%% discretization
M     = 100; % dimension of the model = number of interior grid points
options.model.function     = 'diffusion_filter';
options.model.constants(1) = alpha; % diffusion coefficient
options.model.constants(2) = L; % domain length

x       = linspace(0,L,M+2)';  % grid points including boundary
x_in    = x(2:end-1);
dx      = L/(M+1); % grid spacing, this equals diff(x)

diagonal = ones(M,1);
D        = alpha*(1/(dx^2))*spdiags([diagonal -2*diagonal diagonal],[-1 0 1],M,M);
% time-dependent boundary condition
yD       = zeros(M,1);
yD(1)    = alpha*(1/(dx^2));
yD(end)  = alpha*(1/(dx^2));

options.model.D = D;
options.model.yD = yD;
options.model.x_in = x_in;

%% initial condition
u_start = sin(2*pi*x_in) + sin(8*pi*x_in); % initial condition at interior points

t_start = 0;
t_end   = 1;

N_list  = 10;

% note: dt = (t_end - t_start)/N_list

% settings for solving nonlinear equation
eps_nonlinear = 1e-11;
it_max  = 10;
jacobian = 2; % 1: FD, 2: AD
newton  = 0; % 1: approximate Newton, Jacobian determined once per time step; 2: full Newton
% would be nice to make a special case in the code for linear problems, where the Jacobian only has to be determined once (newton=0)