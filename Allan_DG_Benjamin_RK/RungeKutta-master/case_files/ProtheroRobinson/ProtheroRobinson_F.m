function F = ProtheroRobinson_F(U,t,options,~)
%% Prothero-Robinson problem
% Prothero & Robinson,﻿On the stability and accuracy of one-step methods for solving stiff systems of ordinary differential equations
% the original problem is
% du/dt = g'(t) + lambda*(u - g(t)) 
% exact solution  u = g(t)

% here we take a slightly different problem:
% du/dt = -lambda*(u - cos(t))
% exact solution u(t) = (lambda^2/(lambda^2+1))*(-exp(-lambda*t) + cos(t) + sin(t)/lambda)

        lambda = options.model.constants(1);
        F      = -lambda*(U - cos(t));     

end