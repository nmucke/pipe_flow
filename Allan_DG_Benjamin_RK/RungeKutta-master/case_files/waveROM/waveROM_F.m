function F = waveROM_F(U,~,options)
% d2u/dt2 = nu*d2u/dx^2
% discretized as
% du/dt = v
% dv/dt = D * u
% with periodic boundary conditions

% U = [q;p] = [u; u_t]

ROM_basis = options.model.ROM_basis;

if (strcmp(ROM_basis,'symplectic'))
    % reduced basis:
    % J = options.model.J;
    A = options.model.A;
    
    W = A*U;
    
    
    % approximate state; note that the input U to the function
    % represents the reduced order model solution,
    % whereas W is the approximation to the full order model
    
    J = options.model.J;
    
    D = options.model.D;
    M = 2*options.model.Morig; % size of original system
    
    % Hamiltonian form
    % Z  = zeros(2,2);
    % E  = eye(2);
    % A  = [Z E; -E Z];
    % [dq/dt; dp/dt] = A*[dH/dq; dH/dp];
    % in this case, we have
    % dH/dq = -D*q
    % dH/dp = p
    % this can be simplified to
    % [dq/dt; dp/dt] = [0 I; D 0]*[q;p]; (K*x in paper)
    %                = [0 I; -I 0] * [-D 0; 0 I] * [q;p] (J*L*x in paper)
    %                = [0 I; -I 0] * [-D*q;p] (original Hamiltonian form)% here q=u_t, p=u_tt
    % Z  = spalloc(M/2,M/2,0);
    % E  = speye(M/2);
    % J  = [Z E; -E Z];
    dH = A'*[-D*W(1:M/2); W(M/2+1:end)];
    F  = J*dH; % J2k * A' * L * A * z
    
    % alternative, which uses the original structure F = A'*H(A*U)
    % F  = A'*Jfull*[-D*W(1:M/2); W(M/2+1:end)];
    
    % F = [W(M/2+1:end); D*W(1:M/2)];
    
elseif (strcmp(ROM_basis,'POD'))
    
    V = options.model.V;
    W = V*U;    
    M = options.model.Morig; % number of nodes
    D = options.model.D;    
    F = V'*[W(M+1:end); D*W(1:M)];

else
    
    error('ROM_basis not available');
end

end