function F = Brusselator_F(U,~,options)
% du/dt = A + u^2*v - (B+1)*u + alpha*d2u/dx2
% dv/dt = B*u - u^2*v + alpha*d2v/dx2

% discretized as
% du/dt = A + u^2*v - (B+1)*u + alpha*D*u + uBC
% dv/dt = B*u - u^2*v + alpha*D*v + vBC

% indices of u and v components
ind_u = options.model.ind_u;
ind_v = options.model.ind_v;

Mu = length(ind_u);

A = options.model.constants(3);
B = options.model.constants(4);
D = options.model.D;
uBC = options.model.uBC;
vBC = options.model.vBC;

F =  [A*ones(Mu,1) + U(ind_u).^2 .* U(ind_v) - (B+1)*U(ind_u) + D*U(ind_u) + uBC; ...
      B*U(ind_u) - U(ind_u).^2 .* U(ind_v) + D*U(ind_v) + vBC];

end