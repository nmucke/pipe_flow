function F = Lorenz63_F(U,~,options)

       
        sigma  = options.model.constants(1);
        rho    = options.model.constants(2);
        beta   = options.model.constants(3);

%         F(1,1) = sigma*(U(2)-U(1));
%         F(2,1) = U(1).*(rho-U(3))-U(2);
%         F(3,1) = U(1).*U(2) - beta*U(3);       
        
        F = [sigma*(U(2)-U(1)); ...
             U(1).*(rho-U(3))-U(2); ...
             U(1).*U(2) - beta*U(3)];     

end