function F = simple_F(U,~,options)
% linear ODE
% du/dt = lambda*u

        lambda = options.model.constants(1);
        F      = lambda*U;     

end