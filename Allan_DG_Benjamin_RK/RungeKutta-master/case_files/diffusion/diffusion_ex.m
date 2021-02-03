function u = diffusion_ex(t,options)
% du/dt = nu*d2u/dx^2

alpha = options.model.constants(1); % diffusion coefficient
L     = options.model.constants(2); % domain length

u0      = options.model.u0;

u       = u0.*exp(-pi^2*alpha*t/(L^2));
% + (t.^(p+1) - t_start.^(p+1))/(p+1); 

end