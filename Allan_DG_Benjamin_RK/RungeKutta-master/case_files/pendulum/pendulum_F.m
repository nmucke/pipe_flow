function F = pendulum_F(U,~,options)
% pendulum problem:
% d2 theta/ dt^2 + (g/l) sin(theta) = 0
% first order system: let u = theta, v = dtheta/dt
% du / dt = v
% dv / dt = - (g/l) sin(u)

    l = options.model.constants(1);
    g = options.model.constants(2);
    F = [U(2); -(g/l)*sin(U(1))];

end