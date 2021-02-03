function [rhsu] = AdvecRHS1D(u,time, a)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form field differences at faces
alpha=1;
du = zeros(Nfp*Nfaces,K); 
du(:) = (u(vmapM)-u(vmapP)).*(a*nx(:)-(1-alpha)*abs(a*nx(:)))/2;

% impose boundary condition at x=0
uin = u(vmapO);
uout = u(vmapI);
du (mapI) = (u(vmapI)- uin ).*(a*nx(mapI)-(1-alpha)*abs(a*nx(mapI)))/2;
du (mapO) = (u(vmapO)- uout ).*(a*nx(mapO)-(1-alpha)*abs(a*nx(mapO)))/2;

% compute right hand sides of the semi-discrete PDE
rhsu = -a*rx.*(Dr*u) + LIFT*(Fscale.*(du));
return
