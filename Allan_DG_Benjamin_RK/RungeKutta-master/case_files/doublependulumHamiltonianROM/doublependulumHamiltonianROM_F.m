function F = doublependulumHamiltonianROM_F(U,~,options)
% double pendulum problem:
% U = [U(1); U(2); U(3); U(4)] = [theta1; theta2; p1; p2];

L1 = options.model.constants(1);
L2 = options.model.constants(2);
m1 = options.model.constants(3);
m2 = options.model.constants(4);
g  = options.model.constants(5);

% see e.g. https://diego.assencio.com/?index=e5ac36fcb129ce95a61f8e8ce0572dbf
% or https://mse.redwoods.edu/darnold/math55/DEProj/Sp00/FranScott/finalpaper.pdf
% or http://scienceworld.wolfram.com/physics/DoublePendulum.html
% or http://iopscience.iop.org/article/10.1088/1742-6596/739/1/012066/pdf

% Hamiltonian can be differentiated for example with symbolic toolbox
% syms g m1 m2 L1 L2 p1 p2 theta1 theta2;
% H = (m2*(L2^2)*(p1^2) + (m1+m2)*(L1^2)*(p2^2) - 2*m2*L1*L2*p1*p2*cos(theta1-theta2))/...
%     (2*(L1^2)*(L2^2)*m2*(m1+m2*sin(theta1-theta2).^2)) - ...
%     (m1+m2)*g*L1*cos(theta1) - m2*g*L2*cos(theta2);
% diff(H,p1) % dq1/dt
% diff(H,p2) % dq2/dt
% diff(H,theta1) % -dp1/dt
% diff(H,theta2) % -dp2/dt

% reduced basis:
J = options.model.J;
A = options.model.A;

% approximate state; note that the input U to the function
% represents the reduced order model solution, 
% whereas W is the approximation to the full order model

W = A*U;

theta1 = W(1); %=q1
theta2 = W(2); %=q2
p1     = W(3);
p2     = W(4);

dHdp1 = (2*m2*p1*L2^2 - 2*L1*m2*p2*cos(theta1 - theta2)*L2)/(2*L1^2*L2^2*m2*(m1 + m2*sin(theta1 - theta2)^2));
dHdp2 = (2*p2*(m1 + m2)*L1^2 - 2*L2*m2*p1*cos(theta1 - theta2)*L1)/(2*L1^2*L2^2*m2*(m1 + m2*sin(theta1 - theta2)^2));
dHdq1 = L1*g*sin(theta1)*(m1 + m2) + (p1*p2*sin(theta1 - theta2))/(L1*L2*(m1 + m2*sin(theta1 - theta2)^2)) - (cos(theta1 - theta2)*sin(theta1 - theta2)*((m1 + m2)*L1^2*p2^2 - 2*m2*cos(theta1 - theta2)*L1*L2*p1*p2 + m2*L2^2*p1^2))/(L1^2*L2^2*(m2*sin(theta1 - theta2)^2 + m1)^2);
dHdq2 = L2*g*m2*sin(theta2) - (p1*p2*sin(theta1 - theta2))/(L1*L2*(m1 + m2*sin(theta1 - theta2)^2)) + (cos(theta1 - theta2)*sin(theta1 - theta2)*((m1 + m2)*L1^2*p2^2 - 2*m2*cos(theta1 - theta2)*L1*L2*p1*p2 + m2*L2^2*p1^2))/(L1^2*L2^2*(m2*sin(theta1 - theta2)^2 + m1)^2);

% dq1/dt = dtheta1/dt = (L2*U(3) - L1*U(4)*cos(U(1)-U(2)))/(L1*h3)
% dq2/dt = dtheta2/dt =( -m2*L2*U(3)*cos(U(1)-U(2)) + (m1+m2)*L1*U(4) )/(m2*L2*h3)
% dp1/dt = -(m1+m2)*g*l1*sin(U(1)) - h1 + h2
% dp2/dt = -m2*g*l2*sin(U(2)) + h1 - h2
% h1 = (U(3)*U(4)*sin(U(1)-U(2)))/h3
% h2 =
% sin(2*(U(1)-U(2)))*(m2*(L2^2)*(U(3)^2)+(m1+m2)*(L1^2)*(U(4)^2)-2*m2*L1*L2*U(3)*U(4)*cos(U(1)-U(2)))/
% (2*L1*L2*h3) % there is a factor 2  missing in the second reference?
% h3 = L1*L2*(m1+m2*sin(U(1)-U(2))^2)

% h3 = L1*L2*(m1+m2*sin(U(1)-U(2)).^2);
% h1 = (U(3).*U(4).*sin(U(1)-U(2)))/h3;
% h2 = sin(2*(U(1)-U(2))) .* ...
%     (m2*(L2^2)*(U(3).^2) + (m1+m2)*(L1^2)*(U(4).^2) - 2*m2*L1*L2*U(3).*U(4).*cos(U(1)-U(2))) ./ ...
%     (2*h3.^2); % there is a factor 2 missing in the second and third reference?
% Z  = zeros(2,2);
% E  = eye(2);
% A  = [Z E; -E Z];
dH = [dHdq1; dHdq2; dHdp1; dHdp2];
% F  = A*dH;
F  = J*dH;

% F = [(L2*U(3) - L1*U(4).*cos(U(1)-U(2)))./(L1*h3); ...
%      (-m2*L2*U(3).*cos(U(1)-U(2)) + (m1+m2)*L1*U(4))./(m2*L2*h3); ...
%      -(m1+m2)*g*L1*sin(U(1)) - h1 + h2; ...
%      -m2*g*L2*sin(U(2)) + h1 - h2];

end