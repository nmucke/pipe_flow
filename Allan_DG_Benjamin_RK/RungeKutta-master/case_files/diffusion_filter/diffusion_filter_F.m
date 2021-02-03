function F = diffusion_filter_F(U,t,options)
% du/dt = nu*d2u/dx^2 + f(x,t)
% discretized as
% du/dt = D * u + yD + f
% with zero Dirichlet boundary conditions at x=0 and x=L

% 
x_in = options.model.x_in;
f  = 1 + 4*pi^2*sin(2*pi*x_in) + 64*pi^2*sin(8*pi*x_in);
D  = options.model.D;
% time dependent BC (same for x=0 and x=1)
yD = t*options.model.yD;

F  = D*U + yD + f;
 
end