%% Lotka-Volterra equations

% dimension of problem 
M       = 2; 
options.model.function = 'LotkaVolterra';

options.model.constants(1) = 2;
options.model.constants(2) = 1;
options.model.constants(3) = 1;
options.model.constants(4) = 1;

u_start = [1;1]; % initial condition

N_list  = 100;
t_start = 0;
t_end   = 25;

% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2;