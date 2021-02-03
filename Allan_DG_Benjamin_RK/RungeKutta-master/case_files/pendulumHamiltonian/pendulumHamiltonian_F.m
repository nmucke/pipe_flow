function F = pendulumHamiltonian_F(U,~,options)
% pendulum problem in Hamiltonian formulation
% see e.g. https://cds.cern.ch/record/399399/files/p1.pdf
% H = p^2/(2*l^2) + g*l*(1-cos(theta))
% where p = l^2*dtheta/dt, q=theta
% dq/dt = dtheta/dt = dH/dp = p/l^2
% dp/dt = -dH/dq = -g*l*sin(theta)
% 

    l = options.model.constants(1);
    g = options.model.constants(2);
    F = [U(2)/(l^2); -g*l*sin(U(1))];

end