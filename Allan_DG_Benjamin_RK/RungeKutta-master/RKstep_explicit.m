function unew = RKstep_explicit(un,tn,dt,options)

s = options.RK_method.stages;
M = length(un);

A_RK = options.RK_method.A_RK;
b_RK = options.RK_method.b_RK;
c_RK = options.RK_method.c_RK;

% right hand side evaluations
% F_RK = [u11 u12 ... u1M;
%         u21 ...        ;
%         ...
%         us1 ... ... usM]
F_RK = zeros(s,M);

for kk=1:s

    % intermediate stage value of stage i
    uj   = un + dt*(A_RK(kk,:)*F_RK)'; % use transpose to make result a column vector, size M x 1

    % time level of this stage
    tj   = tn + c_RK(kk)*dt;

    % flux evaluation
    F_RK(kk,:) = getRHS(uj, tj, options);

end

% update solution with the b-coefficients
unew = un + dt*(b_RK'*F_RK)';