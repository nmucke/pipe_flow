function F = wave_F(U,~,options)
% d2u/dt2 = nu*d2u/dx^2
% discretized as
% du/dt = v
% dv/dt = D * u
% with periodic boundary conditions 

% U = [q;p] = [u; u_t]

% 
D = options.model.D;
M = options.model.M;
% Hamiltonian form
%
% [dq/dt; dp/dt] = [0 I; D 0]*[q;p];
% here q=u_t, p=u_tt

F = [U(M/2+1:end); D*U(1:M/2)];

% alternative is the Hamiltonian form:
% Z  = zeros(2,2);
% E  = eye(2);
% A  = [Z E; -E Z];
% [dq/dt; dp/dt] = A*[dH/dq; dH/dp];
% in this case, we have
% dH/dq = -D*q
% dH/dp = p
% this can be written as
% [dq/dt; dp/dt] = [0 I; D 0]*[q;p]; (K*x in paper)
%                = [0 I; -I 0] * [-D 0; 0 I] * [q;p] (J*L*x in paper)
%                = [0 I; -I 0] * [-D*q;p] (original Hamiltonian form)
% here q=u_t, p=u_tt
% Z  = spalloc(M/2,M/2,0);
% E  = speye(M/2);
% J  = [Z E; -E Z];
% dH = [-D*U(1:M/2); U(M/2+1:end)];
% F  = J*dH; 
 
end