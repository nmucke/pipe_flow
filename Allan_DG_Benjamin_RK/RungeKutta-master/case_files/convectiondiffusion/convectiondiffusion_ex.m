function u = convectiondiffusion_ex(t,options)
% du/dt = -c*du/dx + alpha*d2u/dx^2
% discretized as
% du/dt = A*u = -C*u + D*u
% with periodic boundary conditions 
% and IC sin(2*pi*x/L)
% exact solution is sin(2*pi*(x-c*t)/L)*exp(-pi^2*alpha*t/(L^2))

alpha = options.model.constants(1); % diffusion coefficient
% L     = options.model.constants(2); % domain length
c     = options.model.constants(3); % convection coefficient
x_in  = options.model.x_in;
% u0    = options.model.u0;

u      = sin(pi*(x_in-c*t)).*exp(-pi^2*alpha*t);


end