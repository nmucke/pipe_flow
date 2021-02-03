function rhs = NSWERHS_implicit(un,time)
%
% function [rhsh,rhshu] = NSWERHS1D(h,hu,time,g)
% Purpose : Evaluate RHS flux in 1D NSWE
%
% By Allan P. Engsig-Karup, apek@imm.dtu.dk.
Globals1D;

g = 9.81;

h = un(1:end/2);
hu = un(end/2+1:end);

h = reshape(h,N+1,K);
hu = reshape(hu,N+1,K);


% The velocity
u = hu./h;

% The flux function
F1 = hu;
F2 = h.*u.^2+1/2*g*h.^2;

% The eigenvalues of the flux jacobian matrix
lambda = max(abs([u(:)+sqrt(g*h(:)/2) ; u(:)-sqrt(g*h(:)/2)]));

% The factor for the flux jump
C = max(lambda);

% The numerical LF flux
fstar1 = zeros(Nfp*Nfaces,K); 
fstar2 = zeros(Nfp*Nfaces,K); 
fstar1(:) = 1/2*(F1(vmapM)+F1(vmapP))+C/2*nx(:).*(h(vmapM)-h(vmapP)); 
fstar2(:) = 1/2*(F2(vmapM)+F2(vmapP))+C/2*nx(:).*(hu(vmapM)-hu(vmapP));

% The total flux
df1 = zeros(Nfp*Nfaces,K); 
df2 = zeros(Nfp*Nfaces,K); 
df1(:) = nx(:).*(F1(vmapM)-fstar1(:)); 
df2(:) = nx(:).*(F2(vmapM)-fstar2(:)); 

% Compute right hand sides of the semi?discrete PDE 
rhsh  = -rx.*(Dr*(F1)) + LIFT*(Fscale.*df1);
rhshu = -rx.*(Dr*(F2)) + LIFT*(Fscale.*df2); 

rhs = [rhsh rhshu];
rhs = rhs(:);
return

