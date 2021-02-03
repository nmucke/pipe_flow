function F = LotkaVolterra_F(U,~,options)
% Lotka-Volterra predator-prey

alpha     = options.model.constants(1);
beta      = options.model.constants(2);
gamma     = options.model.constants(3);
delta     = options.model.constants(4);

% see for example Butcher:
% y' = z
% z' = mu*(1-y^2)*y' - y = mu*(1-y^2)*z - y
F = [U(1).*(alpha - beta*U(2)); ...
     U(2).*(delta*U(1) - gamma)];

end
