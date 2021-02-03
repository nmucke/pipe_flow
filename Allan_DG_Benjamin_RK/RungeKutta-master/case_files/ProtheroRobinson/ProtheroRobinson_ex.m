function u = ProtheroRobinson_ex(t,options)
% Prothero & Robinson,﻿On the stability and accuracy of one-step methods for solving stiff systems of ordinary differential equations
% the original problem is
% du/dt = g'(t) + lambda*(u - g(t)) 
% exact solution  u = g(t)

% here we take a slightly different problem:
% du/dt = -lambda*(u - cos(t))
% exact solution u(t) = (lambda^2/(lambda^2+1))*(-exp(-lambda*t) + cos(t) + sin(t)/lambda)

lambda  = options.model.constants(1);
% u0      = options.model.u0;
% t_start = options.time.t_start;

u       = (lambda^2/(lambda^2+1))*(-exp(-lambda*t) + cos(t) + sin(t)/lambda); 

end