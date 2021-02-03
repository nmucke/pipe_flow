%% Chemical reaction equations

% dimension of problem 
M       = 3; 
options.model.function = 'reaction';

options.model.constants(1) = 0.04;
options.model.constants(2) = 1e4;
options.model.constants(3) = 3e7;

u_start = [1;0;0]; % initial condition

N_list  = 100;
t_start = 0;
t_end   = 1;

% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2;