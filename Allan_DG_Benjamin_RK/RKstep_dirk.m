function [unew, nsolves, res] = RKstep_dirk(un,tn,dt)
% one time step with diagonally implicit Runge-Kutta (DIRK) method
% starting from un to unew starting from tn with time step dt
% nsolves is number of solves of MxM system, with M the dimension of un

M    = length(un);

Globals1D

eps_nonlinear = 1e-5;
it_max = 50;

%A_RK = options.RK_method.A_RK;
%b_RK = options.RK_method.b_RK;
%c_RK = options.RK_method.c_RK;

s=3; r=1.7588;
g=0.5*(1-cos(pi/18)/sqrt(3)-sin(pi/18));
q=(0.5-g)^2;
A_RK=[g     0    0
    0.5-g g    0
    2*g  1-4*g g];
b_RK=[1/(24*q) 1-1/(12*q) 1/(24*q)]';
c_RK=sum(A_RK,2);

% 
%Fname = [NSWE1D_implicit];

I_M  = speye(M);

% right hand side evaluations
% F_RK = [u11 u12 ... u1M;
%         u21 ...        ;
%         ...
%         us1 ... ... usM]
F_RK = zeros(s,M);

nsolves = 0;

for kk=1:s
    
    % starting guess for intermediate stages
    uj    = un;
    % time level of this stage
    tj    = tn + c_RK(kk)*dt;
    % initialize right-hand side of the current stage
    F_RK(kk,:) = NSWERHS_implicit(uj,tj);
    %F_RK(kk,:) = pipe_flow_rhs_implicit(uj,tj);

    % initialize residual of current stage
    res        = - (uj - un)/dt + (A_RK(kk,:)*F_RK)';
    
    % iteration counter
    k     = 0;
    

    %Jk = Jacobian_single(@(u,t) pipe_flow_rhs_implicit(u,t),un,tn);
    Jk = Jacobian_single(@(u,t) NSWERHS_implicit(u,t),un,tn);


    Q = I_M/dt-A_RK(kk,kk)*Jk;
    % determine LU decomposition
    [L,U] = lu(Q);
    
    while (max(abs(res))>eps_nonlinear)
        
        dUj = U\(L\res);
        
        % update solution vector
        uj  = uj + dUj;
        % update iteration counter
        k   = k+1;
        
        % evaluate rhs for next iteration and check residual based on
        % computed Uj
        %F_RK(kk,:) = pipe_flow_rhs_implicit(uj,tj);
        F_RK(kk,:) = NSWERHS_implicit(uj,tj);
        res        = - (uj - un)/dt + (A_RK(kk,:)*F_RK)';
        
        % number of iterations nonlinear solver
        nsolves    = nsolves+1;
        
        if (k>it_max)
            error(['Newton not converged in ' num2str(it_max) ' iterations']);
        end
        
    end
    
    
    h = SlopeLimitN(reshape(uj(1:end/2),N+1,K));
    hu = SlopeLimitN(reshape(uj(end/2+1:end),N+1,K));
    
    uj = [h,hu];
    uj = uj(:);
    
    %F_RK(kk,:) = pipe_flow_rhs_implicit(uj,tj);
    F_RK(kk,:) = NSWERHS_implicit(uj,tj);
        
    
end

% update solution with the b-coefficients
unew = un + dt*(b_RK'*F_RK)';

h = SlopeLimitN(reshape(unew(1:end/2),N+1,K));
hu = SlopeLimitN(reshape(unew(end/2+1:end),N+1,K));

unew = [h,hu];
unew = unew(:);

% residual
res = max(abs(res));

end