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

%% discretization
M     = 4; % dimension of the model = number of interior grid points
K     = 200; 
options.model.function     = 'NSWE';

x_start = 0;
x_end   = L;

% transfer to structure which will be available in functions
options.model.x = x_in;
options.model.x_start = x_start;
options.model.x_end = x_end;
options.model.M = M;
options.model.K = K;


%% initial condition
t_start = 0;
t_end   = 1;

u_start = 

% number of time steps
N_list  = [1000];% 20 40 80];

% note: dt = (t_end - t_start)/N_list

% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2; % 1: FD, 2: AD
newton  = 2;