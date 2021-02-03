function [operators] = GlobalDGOperators(BCtype)
Globals1D;

Globals1D;
% Compute linear operators using global strong form operators
% mass matrix
invM = V*V';
M = inv(invM);
invMk = 2/(x(end ,1)-x(1,1))*invM;
invM_G = kron(eye(K),invMk);
% differentiation matrix
Dx = Dr;
% stiffness matrix
S = M*Dx;
% flux elemental operators
E = zeros(Np); F = zeros(Np); G = zeros(Np); H = zeros(Np);
E(1,1) = 1; F(1,Np) = 1; G(Np ,1) = 1; H(Np ,Np) = 1;
supdiag = spdiags(ones(K,1) ,1,K,K);
subdiag = spdiags(ones(K,1) ,-1,K,K);
% global average operator - weak form
E_G = kron(subdiag , -0.5*F) + kron(speye(K) ,0.5*(H-E)-S') + kron( - ...
supdiag ,0.5*G);
% global average operator - strong form
G_G = kron(subdiag , -0.5*F) + kron(speye(K),S+0.5*(E-H)) + kron( - ...
supdiag ,0.5*G);
% global jump operator - Lax -Friedrich
F_G = kron(subdiag ,-F) + kron(speye(K),E+H) + kron(supdiag ,-G);

% global average operator - weak form
E_G (1:Np ,end -Np+1: end) = -0.5*F;
E_G(end -Np+1:end ,1:Np) = 0.5*G;
% global average operator - strong form
G_G (1:Np ,end -Np+1: end) = -0.5*F;
G_G(end -Np+1:end ,1:Np) = 0.5*G;
% global jump operator - Lax -Friedrich
F_G (1:Np ,end -Np+1: end) = -F;
F_G(end -Np+1:end ,1:Np) = -G;

% save operators in output structure
operators.invM = invM_G;
operators.E = E_G;
operators.F = F_G;
operators.G = G_G;
end