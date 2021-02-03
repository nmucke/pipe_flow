function F = doublependulum_F(U,~,options)
% double pendulum problem:
% U = [U(1); U(2); U(3); U(4)] = [theta1; theta2; p1; p2;

L1 = options.model.constants(1);
L2 = options.model.constants(2);
m1 = options.model.constants(3);
m2 = options.model.constants(4);
g  = options.model.constants(5);

% see e.g. https://diego.assencio.com/?index=e5ac36fcb129ce95a61f8e8ce0572dbf
% or https://mse.redwoods.edu/darnold/math55/DEProj/Sp00/FranScott/finalpaper.pdf
% or http://scienceworld.wolfram.com/physics/DoublePendulum.html
% or http://iopscience.iop.org/article/10.1088/1742-6596/739/1/012066/pdf
% dtheta1/dt = (L2*U(3) - L1*U(4)*cos(U(1)-U(2)))/(L1*h3)
% dtheta2/dt =( -m2*L2*U(3)*cos(U(1)-U(2)) + (m1+m2)*L1*U(4) )/(m2*L2*h3)
% dp1/dt = -(m1+m2)*g*l1*sin(U(1)) - h1 + h2
% dp2/dt = -m2*g*l2*sin(U(2)) + h1 - h2
% h1 = (U(3)*U(4)*sin(U(1)-U(2)))/h3
% h2 =
% sin(2*(U(1)-U(2)))*(m2*(L2^2)*(U(3)^2)+(m1+m2)*(L1^2)*(U(4)^2)-2*m2*L1*L2*U(3)*U(4)*cos(U(1)-U(2)))/
% (2*L1*L2*h3) % there is a factor 2  missing in the second reference?
% h3 = L1*L2*(m1+m2*sin(U(1)-U(2))^2)

h3 = L1*L2*(m1+m2*sin(U(1)-U(2)).^2);
h1 = (U(3).*U(4).*sin(U(1)-U(2)))/h3;
h2 = sin(2*(U(1)-U(2))) .* ...
    (m2*(L2^2)*(U(3).^2) + (m1+m2)*(L1^2)*(U(4).^2) - 2*m2*L1*L2*U(3).*U(4).*cos(U(1)-U(2))) ./ ...
    (2*h3.^2); % there is a factor 2 missing in the second and third reference?

F = [(L2*U(3) - L1*U(4).*cos(U(1)-U(2)))./(L1*h3); ...
     (-m2*L2*U(3).*cos(U(1)-U(2)) + (m1+m2)*L1*U(4))./(m2*L2*h3); ...
     -(m1+m2)*g*L1*sin(U(1)) - h1 + h2; ...
     -m2*g*L2*sin(U(2)) + h1 - h2];

end