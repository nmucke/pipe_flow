function F = reaction_F(U,~,options)
% REACTION_F
% Chemical reaction equations

alpha     = options.model.constants(1);
beta      = options.model.constants(2);
gamma     = options.model.constants(3);

% see for example Hairer:

F = [-alpha*U(1) + beta*U(2)*U(3) ; ...
      alpha*U(1) - beta*U(2)*U(3) - gamma*U(2).^2; ...
      gamma*U(2).^2];

end
