function rhs = pipe_flow_rhs_implicit(un,time)
%
% function [rhsh,rhshu] = NSWERHS1D(h,hu,time,g)
% Purpose : Evaluate RHS flux in 1D NSWE
%
% By Allan P. Engsig-Karup, apek@imm.dtu.dk.
Globals1D;
global c p0 rho0

h = un(1:end/2);
hu = un(end/2+1:end);

h = reshape(h,N+1,K);
hu = reshape(hu,N+1,K);

g = 9.81;

% The velocity
u = hu./h;

lam = max(abs([u + c, u - c]));

C = max(lam);

% The flux function

pressure = c*c*(h-rho0) + p0;

F1 = hu;
F2 = h.*u.^2 + pressure;

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

% uIn = 0.;
% uOut = 0.;
% huIn = uIn.*h(vmapI);
% huOut = uOut.*h(vmapO);
% 
% hFluxIn = huIn;
% huFluxIn = h(vmapI) .* uIn.^2 + pressure(vmapI);
% hFluxOut = huIn;
% huFluxOut = h(vmapO) .* uOut.^2 + pressure(vmapO);
% 
% 
% df1(mapI) = nx(mapI)*(F1(vmapI)-hFluxIn)/2;
% df2(mapI) = nx(mapI)*(F2(vmapI)-huFluxIn)/2 - C*(hu(vmapI)-huIn);
% 
% df1(mapO) = nx(mapO)*(F1(vmapO)-hFluxOut)/2;
% df2(mapO) = nx(mapO)*(F2(vmapO)-huFluxOut)/2 - C*(hu(vmapO)-huOut);

% Compute right hand sides of the semi?discrete PDE 


rhsh  = -rx.*(Dr*(F1)) + LIFT*(Fscale.*df1);
rhshu = -rx.*(Dr*(F2)) + LIFT*(Fscale.*df2); 

rhs = [rhsh rhshu];
rhs = rhs(:);
return




















