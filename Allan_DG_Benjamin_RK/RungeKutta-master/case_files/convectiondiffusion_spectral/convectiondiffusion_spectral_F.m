function F = convectiondiffusion_spectral_F(U,~,options)
% du/dt = -c*du/dx + alpha*d2u/dx^2
% discretized as
% du/dt = A*u = -C*u + D*u
% with periodic boundary conditions 
 
A = options.model.A;
F = A*U;
 
end