%% polynomial in time
% du/dt = exp(-lambda*t^2)
% exact solution:
% erf(t) = 2/sqrt(pi) int  exp(-t^2) dt
% int exp(-lambda t^2) dt = int exp(- sqrt(lambda)^2 t^2) dt = 
% int exp( - (sqrt(lambda)*t)^2) dt = int exp( - (sqrt(lambda)*t)^2)
% sqrt(lambda)*dt/sqrt(lambda) = 1/sqrt(lambda) exp ( - tau^2) dtau; (note tau
% = sqrt(lambda)*t )
% = (1/2)*sqrt(pi/lambda)*erf(lambda*t)

% useful for checking order of accuracy since methods of order s should
% integrate polynomials up to order s-1 exactly

M                          = 1; % dimension of the model
options.model.function     = 'Trefethen1';
options.model.constants(1) = 100; % lambda
u_start                    = 0; % initial condition

t_start = -0.5;
t_end   = 1;

N_list  = [1 2 4 8 16 32 64 128];
% N_list = [100 200 400]*10;

% note: dt = (t_end - t_start)/N_list

% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian  = 2;