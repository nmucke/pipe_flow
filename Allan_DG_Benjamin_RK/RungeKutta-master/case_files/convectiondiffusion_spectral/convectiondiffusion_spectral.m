%% PDE:
% du/dt = -c*du/dx + alpha*d2u/dx^2
% discretized as
% du/dt = A*u = -C*u + D*u
% with a Fourier spectral method and periodic boundary conditions 

addpath('/Users/sanderse/Dropbox/work/Programming/libs/chebfun-master/');

%% physical values
% L     = 1; % length of domain
alpha = 0.1; % diffusion coefficient
c     = 1; % convection coefficient


%% discretization
M     = 25; % dimension of the model = number of interior grid points
options.model.function     = 'convectiondiffusion_spectral';
options.model.constants(1) = alpha; % diffusion coefficient
% options.model.constants(2) = L; % domain length
options.model.constants(3) = c; % convection coefficient

x       = linspace(-1,1,M+1)';  % grid points including boundary
x_in    = x(1:end-1); % periodic bc, leave out last point
dx      = (x(end)-x(1))/M; % grid spacing, this equals diff(x)

options.model.x_in= x_in;

% use Chebfun for building discretization:
X = [-1,1];

% diffusion:
D = alpha*diffmat(M,2,'periodic',X); % MxM Fourier differentiation matrix of order 2 on interval X

% convection:
C = c*diffmat(M,1,'periodic',X); % MxM Fourier differentiation matrix of order 1 on interval X

options.model.A = -C + D;

%% initial condition
u_start = sin(pi*x_in); % initial condition

t_start = 0;
t_end   = 2.5;

N_list  = 200;

% note: dt = (t_end - t_start)/N_list

% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2; % 1: FD, 2: AD
newton  = 1; % would be nice to make a special case in the code for linear problems, where the Jacobian only has to be determined once (newton=0)