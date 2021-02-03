function [unew, nsolves, res] = RKstep_dirk(un,tn,dt,options)
% one time step with diagonally implicit Runge-Kutta (DIRK) method
% starting from un to unew starting from tn with time step dt
% nsolves is number of solves of MxM system, with M the dimension of un

s    = options.RK_method.stages;
M    = length(un);

eps_nonlinear = options.settings.eps_nonlinear;
it_max = options.settings.it_max;

A_RK = options.RK_method.A_RK;
b_RK = options.RK_method.b_RK;
c_RK = options.RK_method.c_RK;

% 
Fname = [options.model.function '_F'];

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
    F_RK(kk,:) = getRHS(uj,tj,options);
    % initialize residual of current stage
    res        = - (uj - un)/dt + (A_RK(kk,:)*F_RK)';
    
    % iteration counter
    k     = 0;
    
    if (options.settings.newton == 0 || options.settings.newton == 1)
        if (options.settings.newton == 0)
            Jk = options.settings.J;
        elseif (options.settings.newton == 1) % approximate Newton
            % Jacobian based on current solution un
            Jk = Jacobian_single(Fname,un,tn,options);
            % form iteration matrix, which is now fixed during this stage
        end
        Q = I_M/dt-A_RK(kk,kk)*Jk;
        % determine LU decomposition
        [L,U] = lu(Q);
    end
    
    while (max(abs(res))>eps_nonlinear)
        
        if (options.settings.newton == 0 || options.settings.newton ==  1) 
            % approximate Newton
            % re-use the LU decomposition
            dUj = U\(L\res);
            
        elseif (options.settings.newton == 2)
            % full Newton
            J   = Jacobian(Fname,uj,tj,options,1);        
            % form iteration matrix
            Q   = I_M/dt-A_RK(kk,kk)*J;
            % get change
            dUj = Q\res;
        end
        
        % update solution vector
        uj  = uj + dUj;
        % update iteration counter
        k   = k+1;
        
        % evaluate rhs for next iteration and check residual based on
        % computed Uj
        F_RK(kk,:) = getRHS(uj,tj,options);
        res        = - (uj - un)/dt + (A_RK(kk,:)*F_RK)';
        
        % number of iterations nonlinear solver
        nsolves    = nsolves+1;
        
        if (k>it_max)
            error(['Newton not converged in ' num2str(it_max) ' iterations']);
        end
        
    end
    
    
end

% update solution with the b-coefficients
unew = un + dt*(b_RK'*F_RK)';

% residual
res = max(abs(res));