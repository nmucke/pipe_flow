%% PDE:
% du/dt + u du/dx = nu*d2u/dx^2
% or equivalently
% du/dt + d (0.5*u^2) /dx = nu*d2u/dx^2


%% physical values
L     = 1; % length of domain
nu    = 0.001; % diffusion coefficient


%% 
% u = spin('Burg'); %,'plot','off');

%%
dom = [-1 1]; 
x = chebfun('x',dom); 
tspan = 0:10:100;
Ngrid = 100;
dt = 0.01;
S = spinop(dom,tspan);
S.lin = @(u) 0.01*diff(u,2);
S.nonlin = @(u) -.5*diff(u.^2); % spin cannot parse "u.*diff(u)"
% S.init = (1-x^2)*exp(-30*(x+1/2)^2);
% S.init = sin(x);
S.init = 0.001*2*pi*sin(pi*x)./(2+cos(pi*x));
u = spin(S,Ngrid,dt); %,'plot','off');
% waterfall(u)
% plot(S.init)
% hold on
% plot(u{end})
% ylim([-4 4]), hold off
% FS = 'fontsize'; text(42,3.4,'t=0 and t=100',FS,12)
hold on

%% 
% exact solution:
t=tspan(end);
nu = 0.001;
u = nu*2*pi*exp(-(pi^2)*t*nu)*sin(pi*x)./(2+exp(-(pi^2)*t*nu)*cos(pi*x));
plot(x,u)

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