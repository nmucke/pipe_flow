function F = Burgers_F(U,t,options)
% BURGERS_F right-hand side of Burgers equation
% du/dt = F(u(t),t) =  - d/dx(0.5* u^2) + nu*d2u/dx^2
% = -C*(u.^2) + D*u + BC
% see Rang, JCAM 2015, Improved traditional Rosenbrock–Wanner methods for stiff
% ODEs and DAEs


x_start = options.model.x_start;
x_end   = options.model.x_end;

D = options.model.D;
C = options.model.C;
A = options.model.A;

% M = options.model.M;

% use exact solution for the boundary conditions
uBCL     = Burgers_ex(t,options,x_start);
uBCR     = Burgers_ex(t,options,x_end);
% uBCL = 0;
% uBCR = 1;

% we use the discretization matrices to construct the BC
% uBC      = zeros(M,1);
% uBC(1)   = (+C(1,2)*uBCL^2 + D(1,2)*uBCL); % note the plus before C(1,2), which is -*-
% uBC(M)   = (+C(M,M-1)*uBCR^2 + D(M,M-1)*uBCR);

Utot     = [uBCL;U;uBCR];

% right-hand side
F        = -C*((A*Utot).^2) + D*Utot;
 
end