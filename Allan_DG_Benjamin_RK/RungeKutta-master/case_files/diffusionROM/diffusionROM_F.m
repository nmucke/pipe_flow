function F = diffusionROM_F(U,~,options)
% du/dt = nu*d2u/dx^2
% discretized as
% du/dt = D * u
% with zero Dirichlet boundary conditions at x=0 and x=L

% D is the ROM diffusion matrix, which is precomputed
D = options.model.D;
F = D*U;


end