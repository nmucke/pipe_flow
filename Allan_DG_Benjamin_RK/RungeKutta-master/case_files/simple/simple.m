%% linear ODE
% du/dt = lambda*u
% exact solution u(t) = u(t0)*exp(lambda*t);

M        = 1;
lambda   = -1;
u_start  = 1;

options.model.function     = 'simple';
options.model.constants(1) = lambda;

t_start = 0;
t_end   = 1;

N_list  = [10 20 40 80 160];


% settings for solving nonlinear equation
eps_nonlinear = 1e-6;
it_max  = 20;
jacobian = 1;
newton = 1;