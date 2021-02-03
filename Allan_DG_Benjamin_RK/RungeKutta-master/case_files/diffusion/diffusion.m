%% PDE:
% du/dt = nu*d2u/dx^2
% discretized as
% du/dt = D * u
% with zero Dirichlet boundary conditions at x=0 and x=L

%% physical values
L     = 10; % length of domain
alpha = 1;


%% discretization
M     = 200; % dimension of the model = number of interior grid points
options.model.function     = 'diffusion';
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
newton  = 1; % would be nice to make a special case in the code for linear problems, where the Jacobian only has to be determined once (newton=0)