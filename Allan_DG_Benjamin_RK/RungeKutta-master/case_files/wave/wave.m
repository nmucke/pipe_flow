%% PDE:
% d2u/dt2 = nu*d2u/dx^2
% discretized as
% du/dt = v
% dv/dt = D * u
% with zero Dirichlet boundary conditions at x=0 and x=L

%% physical values
L     = 1; % length of domain
c     = 0.1; % wave speed


%% discretization
M     = 500; %  number of interior grid points
options.model.function     = 'wave';
options.model.constants(1) = c; % diffusion coefficient
options.model.constants(2) = L; % domain length

x       = linspace(0,L,M+1)';  % grid points
x_in    = x(1:end-1); % last point equals first point due to periodic BC
dx      = L/(M+1); % grid spacing, this equals diff(x)

diagonal = ones(M,1);
D        = (1/(dx^2))*spdiags([diagonal -2*diagonal diagonal],[-1 0 1],M,M);
% periodic BC
D(1,end) = 1/(dx^2);
D(end,1) = 1/(dx^2);
options.model.D = (c^2)*D;

%% initial condition

% u_start = sin(pi*x_in/L); % initial condition
s       = 10*abs(x_in-0.5);
q_start = (1-(3/2)*(s.^2)+(3/4)*(s.^3)).*(s>=0 & s<=1) + ...
           (1/4)*(2-s).^3.*(s>=1 & s<2);
p_start = zeros(M,1);
u_start = [q_start; ...
           p_start];
              
M = 2*M; % dimension of the model

t_start = 0;
t_end   = 50;

N_list  = t_end/0.01;

% note: dt = (t_end - t_start)/N_list

% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2; % 1: FD, 2: AD
newton  = 2;