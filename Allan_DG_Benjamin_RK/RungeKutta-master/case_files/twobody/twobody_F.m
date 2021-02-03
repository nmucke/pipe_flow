function F = twobody_F(U,~,~)
% See Hairer et al., 'Geometrical Numerical Integration':
% d^2/dt^2 q1 = -q1/(q1^2 + q2^2)^(3/2)
% d^2/dt^2 q2 = -q2/(q1^2 + q2^2)^(3/2)
% write as first order system
% d/dt U1 = d/dt q1 = v1
% d/dt U2 = d/dt q2 = v2
% d/dt U3 = d/dt v1 = -q1/(q1^2 + q2^2)^(3/2)
% d/dt U4 = d/dt v2 = -q2/(q1^2 + q2^2)^(3/2)

% U1 = q1
% U2 = q2
% U3 = v1
% U4 = v2

R =  (U(1).^2  + U(2).^2).^(3/2);
F = [U(3);...
     U(4);...
     -U(1)./R;...
     -U(2)./R];

end
