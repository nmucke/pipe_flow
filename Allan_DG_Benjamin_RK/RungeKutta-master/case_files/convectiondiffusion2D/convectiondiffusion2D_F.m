function F = convectiondiffusion2D_F(U,~,options)
% du/dt = -c*du/dx + alpha*d2u/dx^2
% discretized as
% du/dt = A*u = -C*u + D*u + r
% with boundary conditions given by r
 
A  = options.model.A;
bc = options.model.bc.r;
F  = A*U + bc;
 
end