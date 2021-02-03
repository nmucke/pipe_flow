function F = polynomial_F(U,t,options)
% not a real ODE, only depending on time:
% du/dt = t^p
        p = options.model.constants(1);
        F = t.^p;     

end