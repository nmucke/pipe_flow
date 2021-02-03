function F = Trefethen1_F(U,t,options)
% du/dt = exp(-lambda*t^2)

        lambda = options.model.constants(1);
        F = exp(-lambda*t.^2);     

end