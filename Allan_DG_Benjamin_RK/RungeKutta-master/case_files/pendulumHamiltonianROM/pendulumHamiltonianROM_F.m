function F = pendulumHamiltonianROM_F(U,~,options)
% pendulum problem in Hamiltonian formulation
% see e.g. https://cds.cern.ch/record/399399/files/p1.pdf
% H = p^2/(2*l^2) + g*l*(1-cos(theta))
% where p = l^2*dtheta/dt
% dq/dt = dtheta/dt = dH/dp = p/l^2
% dp/dt = -dH/dq = -g*l*sin(theta)

    % reduced basis:
    J = options.model.J;
    A = options.model.A;
    
    % approximate state; note that the input U to the function
    % represents the reduced order model solution, 
    % whereas W is the approximation to the full order model
    
    W = A*U;
    l = options.model.constants(1);
    g = options.model.constants(2);
    % dH = [dH/dq; dH/dp];
    dH = [g*l*sin(W(1)); W(2)/(l^2)];
    F  = J*dH;

end