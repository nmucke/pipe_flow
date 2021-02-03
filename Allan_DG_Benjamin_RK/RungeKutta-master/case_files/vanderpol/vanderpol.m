%% vanderpol:

% dimension of problem 
M       = 2; 
options.model.function = 'vanderpol';
mu      = 5;
options.model.constants(1) = mu;

u_start = [0;1]; % initial condition

N_list  = 200;
t_start = 0;
t_end   = 50;

% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian  = 2;