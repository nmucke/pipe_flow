%% Prothero-Robinson problem
% Prothero & Robinson,﻿On the stability and accuracy of one-step methods for solving stiff systems of ordinary differential equations
% the original problem is
% du/dt = g'(t) + lambda*(u - g(t)) 
% exact solution  u = g(t)

% here we take a slightly different problem:
% du/dt = -lambda*(u - cos(t))
% exact solution u(t) = (lambda^2/(lambda^2+1))*(-exp(-lambda*t) + cos(t) + sin(t)/lambda)


M                          = 1; % dimension of the model
options.model.function     = 'ProtheroRobinson';
options.model.constants(1) = 2000; % lambda
u_start                    = 0; % initial condition

t_start = 0;
t_end   = 1;

N_list  = [10 20 40 80 160 320 640 1280];
% N_list = [100 200 400]*10;

% note: dt = (t_end - t_start)/N_list

% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2; % 1: FD, 2: AD
newton = 2;
