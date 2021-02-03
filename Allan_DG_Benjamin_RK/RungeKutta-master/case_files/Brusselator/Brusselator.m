%% PDE:
% du/dt = A + u^2*v - (B+1)*u + alpha*d2u/dx2
% dv/dt = B*u - u^2*v + alpha*d2v/dx2

% discretized as
% du/dt = A + u^2*v - (B+1)*u + alpha*D*u + uBC
% dv/dt = B*u - u^2*v + alpha*D*v + vBC

%% physical values
L     = 1; % length of domain
alpha = 1/50;
A     = 1;
B     = 3;

%% discretization
Mi     = 200; % number of interior grid points per equation
options.model.function     = 'Brusselator';

options.model.constants(1) = alpha; % diffusion coefficient
options.model.constants(2) = L; % domain length
options.model.constants(3) = A;
options.model.constants(4) = B;

% grid
x       = linspace(0,L,Mi+2)';  % grid points including boundary
x_in    = x(2:end-1);
dx      = L/(Mi+1); % grid spacing, this equals diff(x)

% diffusion matrix
diagonal = ones(Mi,1);
D        = alpha*(1/(dx^2))*spdiags([diagonal -2*diagonal diagonal],[-1 0 1],Mi,Mi);
options.model.D = D;

% boundary conditions
uBCL     = 1;
uBCR     = 1;
vBCL     = 3;
vBCR     = 3;
uBC      = zeros(Mi,1);
uBC(1)   = D(1,2)*uBCL;
uBC(Mi)  = D(1,2)*uBCR;
vBC      = zeros(Mi,1);
vBC(1)   = D(1,2)*vBCL;
vBC(Mi)  = D(1,2)*vBCR;
options.model.uBC = uBC;
options.model.vBC = vBC;

% index of u and v in system of equations
options.model.ind_u = (1:Mi)';
options.model.ind_v = (Mi+1:2*Mi)';

% dimension of the problem
M = 2*Mi;

%% initial condition

% note: u = [u;v]; initial conditions for u and v
u_start = [1+sin(2*pi*x_in/L); 3*ones(Mi,1)];


%% time interval and time steps
t_start = 0;
t_end   = 10;

% number of time steps
N_list  = 50;

% note: dt = (t_end - t_start)/N_list

% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2; % 1: FD, 2: AD
newton  = 2;