function F = getRHS(U,t,options,nstages)
% U is a column vector with current solution
% t is vector with time instance of stages
% F is model rhs evaluated at (U,t)
global rhs_evals;

if (nargin<=3)
    nstages=1;
end

M = options.model.M;
F = zeros(nstages*M,1);

% evaluate F for each stage
for k=1:nstages
    % ind1 and ind2 are the indices in the global vector
    ind1 = M*(k-1) + 1;
    ind2 = M*k;
    
    rhs_evals = rhs_evals+1;
    
    func         = str2func([options.model.function '_F']);
    F(ind1:ind2) = func(U(ind1:ind2),t(k),options);
    
    
end
end
