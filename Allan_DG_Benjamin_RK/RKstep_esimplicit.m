function [unew, nsolves, res] = RKstep_esimplicit(un,tn,dt,options)
% one time step with implicit Runge-Kutta method
% starting from un to unew starting from tn with time step dt
% nsolves is number of solves of sMxsM system, with M the dimension of un,
% and s the number of stages of the RK method

s    = options.RK_method.stages;
M    = length(un);

eps_nonlinear = options.settings.eps_nonlinear;
it_max = options.settings.it_max;


A_RK = options.RK_method.A_RK;
b_RK = options.RK_method.b_RK;
c_RK = options.RK_method.c_RK;

I_sM      = kron(speye(s-1),speye(M));
A_RK_ext  = kron(A_RK(2:end,:),speye(M));
A_RK_in   = kron(A_RK(2:end,2:end),speye(M));
b_RK_ext  = kron(b_RK',speye(M));

%
Fname = [options.model.function '_F'];


% capital U contains all stages and is ordered as [u1;v1;w1;...;u2;v2;w2;...;us;vs;ws];
Un    = kron(ones(s-1,1),un);
% tj = [t1;t2;...;ts]
tj    = tn + c_RK(2:end)*dt;

%
F = zeros(s*M,1);

% F for first stage
F(1:M) = getRHS(un,tn,options,1);

% iteration counter
k     = 0;

% starting guess for intermediate stages
Uj    = Un;
% initialize right-hand side
F(M+1:s*M)  = getRHS(Uj,tj,options,s-1);
% initialize residual
res   = - (Uj - Un)/dt + A_RK_ext*F;% + A_RK1*F1;


if (options.settings.newton == 0) 
   L = options.settings.Jacobian_L; 
   U = options.settings.Jacobian_U; 
end
if (options.settings.newton == 1) % approximate Newton
    % Jacobian based on current solution un
    Jn = Jacobian_single(Fname,un,tn,options);
    % form iteration matrix, which is now fixed during iterations
    Q = (I_sM/dt-kron(A_RK(2:end,2:end),Jn));
    % determine LU decomposition
    [L,U] = lu(Q);
end

while (max(abs(res))>eps_nonlinear)
    
    if (options.settings.newton == 1)
        % approximate Newton
        % re-use the LU decomposition
        dUj = U\(L\res);
        
    elseif (options.settings.newton == 2)
        % full Newton
        J   = Jacobian(Fname,Uj,tj,options,s-1);
        % form iteration matrix
        Q   = I_sM/dt-A_RK_in*J;
        % get change
        dUj = Q\res;
    end
    
    % update solution vector
    Uj  = Uj + dUj;
    % update iteration counter
    k   = k+1;
    
    % evaluate rhs for next iteration and check residual based on
    % computed Uj
    F(M+1:s*M)  = getRHS(Uj,tj,options,s-1);
    res         = - (Uj - Un)/dt + A_RK_ext*F; % + A_RK1*F1;
    
    if (k>it_max)
        error(['Newton not converged in ' num2str(it_max) ' iterations']);
    end
    
end

% solution at new time step with b-coefficients of RK method
unew    = un + dt*b_RK_ext*F;

% number of iterations nonlinear solver
nsolves = k;

% residual
res = max(abs(res));