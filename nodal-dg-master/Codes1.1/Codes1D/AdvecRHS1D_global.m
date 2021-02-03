function [rhsu] = AdvecRHS1D(u,time, a)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection
Globals1D;
op = GlobalDGOperators('periodic');

% compute right hand sides of the semi-discrete PDE
rhsu = -op.invM*(a*op.G+0.5*a*op.F)*u(:);

rhsu = reshape(rhsu,Np,K);
return
