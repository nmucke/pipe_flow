function F = vanderpol_F(U,~,options)
% vanderpol oscillator

mu     = options.model.constants(1);

% see for example Butcher:
% y' = z
% z' = mu*(1-y^2)*y' - y = mu*(1-y^2)*z - y
F    = [U(2); ...
        mu*(1-U(1).^2).*U(2) - U(1)];

end
