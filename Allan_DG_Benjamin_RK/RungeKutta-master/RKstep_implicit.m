function [unew, nsolves, res] = RKstep_implicit(un,tn,dt,options,previous)
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

I_sM      = kron(speye(s),speye(M));
A_RK_ext  = kron(A_RK,speye(M));
b_RK_ext  = kron(b_RK',speye(M));

%
Fname = [options.model.function '_F'];

% capital U contains all stages and is ordered as [u1;v1;w1;...;u2;v2;w2;...;us;vs;ws];
% initialize U with the solution at tn
Un    = kron(ones(s,1),un);

% tj contains the time instances at all stages, tj = [t1;t2;...;ts]
tj    = tn + c_RK*dt;

% iteration counter
k     = 0;

% starting guess for intermediate stages
Uj    = Un;
% initialize right-hand side
F     = getRHS(Uj,tj,options,s);
% initialize residual
res   = - (Uj - Un)/dt + A_RK_ext*F;

if (options.settings.newton == 0) 
   L = options.settings.Jacobian_L; 
   U = options.settings.Jacobian_U; 
end

if (options.settings.newton == 1) % approximate Newton
    % Jacobian based on current solution un 
    Jn = Jacobian_single(Fname,un,tn,options);
    % form iteration matrix, which is now fixed during iterations
    Q = (I_sM/dt-kron(A_RK,Jn));
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
        J   = Jacobian(Fname,Uj,tj,options,s);
        % form iteration matrix
        Q   = I_sM/dt-A_RK_ext*J;
        % get change
        dUj = Q\res;
    end
    

    % update solution vector
    Uj  = Uj + dUj;
    % update iteration counter
    k   = k+1;
    
    % evaluate rhs for next iteration and check residual based on
    % computed Uj
    F     = getRHS(Uj,tj,options,s);
    res   = - (Uj - Un)/dt + A_RK_ext*F;
%     max(abs(res));
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