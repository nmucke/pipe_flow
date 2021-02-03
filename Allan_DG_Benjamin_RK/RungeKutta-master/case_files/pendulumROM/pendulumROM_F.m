function F = pendulumROM_F(U,~,options)
% pendulum problem:
% d2 theta/ dt^2 + (g/l) sin(theta) = 0
% first order system: let u = theta, v = dtheta/dt
% du / dt = v
% dv / dt = - (g/l) sin(u)

    % reduced basis:
    V = options.model.V;
    
    % approximate state; note that the input U to the function
    % represents the reduced order model solution, 
    % whereas W is the approximation to the full order model
    W = V*U;
    l = options.model.constants(1);
    g = options.model.constants(2);
    
    F = V'*[W(2); -(g/l)*sin(W(1))];

end