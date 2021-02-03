%% polynomial in time
% du/dt = t^p
% exact solution u(t) = u(t0)*t^(p+1)/(p+1);
% useful for checking order of accuracy since methods of order s should
% integrate polynomials up to order s-1 exactly

M                          = 1; % dimension of the model
options.model.function     = 'polynomial';
options.model.constants(1) = 6; % polynomial order
u_start                    = 1; % initial condition

t_start = 0;
t_end   = 1;

N_list  = [10 20 40 80];
% N_list = [100 200 400]*10;

% note: dt = (t_end - t_start)/N_list

% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian  = 2;