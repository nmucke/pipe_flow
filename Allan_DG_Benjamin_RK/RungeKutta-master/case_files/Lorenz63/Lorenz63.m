%% Lorenz 63 nonlinear ODE


t_start = 0;
t_end   = 40;

N_list = [200 ];
% N_list = [100 200 400]*10;
% N_list = (t_end-t_start)/dt;

M       = 3; % dimension of problem 
options.model.function = 'Lorenz63';
options.model.constants(1)  = 10; % sigma
options.model.constants(2)  = 28; % rho
options.model.constants(3)  = 8/3; % beta
u_start = [1;0;0]; % initial condition


% settings for solving nonlinear equation
eps_nonlinear = 1e-6;
it_max  = 10;
jacobian = 2;
newton  = 2;
