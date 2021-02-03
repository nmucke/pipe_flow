function J = Jacobian( Fname, U, t, options, nstages)
%JACOBIAN
% construct entire Jacobian J=dF/dU for all stages of RK method
% Fname is a string with function name
% U is state around which linearization is calculated

% note: this is an expensive implementation, only for testing purposes
if (nargin<=4)
    nstages=1;
end

% size of state vector
N   = length(U);
M   = N/nstages;

% allocate space for Jacobian
% J    = zeros(N,N);
J = [];

% evaluate J: J1=dF/dU(U1); J2=dF/dU(U2); J3=dF/dU(U3)
% put these as diagonal blocks in matrix J
for k=1:nstages
    % ind1 and ind2 are the indices in the global vector
    ind1 = M*(k-1) + 1;
    ind2 = M*k;
    
    % construct Jacobian k
    Jk = Jacobian_single( Fname, U(ind1:ind2), t(k), options);
    J  = blkdiag(J,Jk);
    
    % alternative:
    % J(ind1:ind2,ind1:ind2) = Jacobian_single( Fname, U(ind1:ind2), t(k), options);
    
end

% make the Jacobian sparse
if (~issparse(J))
    J = sparse(J);
end

end

