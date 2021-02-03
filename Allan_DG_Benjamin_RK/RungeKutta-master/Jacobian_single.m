function J = Jacobian_single(Fname, U, t, options)
%JACOBIAN
% construct Jacobian J=dF/dU by finite differences or automatic
% differentiation
% Fname is a string with function name
% U is state around which linearization is calculated

% note: this is an expensive implementation, only for testing purposes
% if (nargin<=4)
%     nstages=1;
% end

% size of state vector
N   = length(U);

% make F callable as a function
% Fname should exist in the path
func = str2func(Fname);

    
% perturbation determined by machine-epsilon
eps_jacobian = sqrt(eps);

% call F
F    = func(U,t,options);

% size of F can be different from size of W (non-square Jacobians)
M    = length(F);

% allocate space for Jacobian
J    = zeros(M,N);

% vector used for perturbation
z    = zeros(N,1);

for col=1:N

    pert        = z;
    pert_jac    = eps_jacobian*max(abs(U(col)),1);
    pert(col)   = pert_jac;

    Upert     = U + pert;

    Fpert     = func(Upert,t,options);

    J(:,col)  = (Fpert - F) / pert_jac;

end


% make the Jacobian sparse
if (~issparse(J))
    J = sparse(J);
end

end

