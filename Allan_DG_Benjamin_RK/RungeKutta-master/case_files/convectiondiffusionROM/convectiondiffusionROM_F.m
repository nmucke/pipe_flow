function F = convectiondiffusionROM_F(U,~,options)
% du/dt = -c*du/dx + alpha*d2u/dx^2
% discretized as
% du/dt = A*u = -C*u + D*u
% with periodic boundary conditions 
 
% A is the ROM diffusion matrix, which is precomputed
A = options.model.A;
F = A*U;
 
end