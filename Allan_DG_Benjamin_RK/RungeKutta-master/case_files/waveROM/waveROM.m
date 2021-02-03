%% PDE:
% d2u/dt2 = c^2*d2u/dx^2
% discretized as
% du/dt = v
% dv/dt = D * u
% with zero Dirichlet boundary conditions at x=0 and x=L

%% physical values
L     = 1; % length of domain
c     = 0.1; % wave speed


%% discretization
M     = 500; %  number of interior grid points of original (high-fidelity) solution
options.model.function     = 'waveROM';
options.model.constants(1) = c; % diffusion coefficient
options.model.constants(2) = L; % domain length

x       = linspace(0,L,M+1)';  % grid points
x_in    = x(1:end-1); % last point equals first point due to periodic BC
dx      = L/(M+1); % grid spacing, this equals diff(x)

diagonal = ones(M,1);
D        = (1/(dx^2))*spdiags([diagonal -2*diagonal diagonal],[-1 0 1],M,M);
% periodic BC
D(1,end) = 1/(dx^2);
D(end,1) = 1/(dx^2);
options.model.D = (c^2)*D;

%% initial condition

s       = 10*abs(x_in-0.5);
q_start = (1-(3/2)*(s.^2)+(3/4)*(s.^3)).*(s>=0 & s<=1) + ...
          ((1/4)*(2-s).^3).*(s>=1 & s<2);
p_start = zeros(M,1);
u_start = [q_start; p_start];


%% time integration settings
t_start = 0;
t_end   = 1;

N_list  = t_end/0.01;

% note: dt = (t_end - t_start)/N_list

% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2; % 1: FD, 2: AD
newton  = 2;


%% construct reduced basis

% skip relates to number of snapshots used; Nsnaps = Noriginal/skip
skip   = 50;
% retain k leading singular values, k<=Ndim
k      = 40;
% load snapshots from high-fidelity sim
snapshots = load('wave_case2.mat');
% choose for symplectic or standard POD
ROM_basis = 'symplectic';

if (strcmp(ROM_basis,'symplectic'))
    % the snapshot matrix is constructed based on u_plot:
    % u_plot = [q(t1) q(t2) ... q(tN);
    %           p(t1) p(t2) ... p(tN)];
    % with size Ndim x Nt, where Ndim is the dimension of the vector of
    % unknowns (= number of ODEs); which is 2 in this case
    % following Peng & Mohsen, we construct the matrix
    % S = Sq + i * Sp
    
    Morig  = M; % original system size
    I      = sqrt(-1);
    indx_q = 1:M;
    indx_p = M+1:2*M;
    S      = snapshots.u_plot(indx_q,1:skip:end) + I*snapshots.u_plot(indx_p,1:skip:end);
    
    % compute the economic SVD of S: S = K*Sigma*L'
    % K is Ndim x Ndim, Sigma is Ndim x Ndim, L is Nt x Ndim
    % in this case Ndim<<Nt
    [K,Sigma,L] = svd(S,0);
    
    % construct reduced basis from K, such that V'*V = I
    V   = K(:,1:k);
    phi = real(V);
    psi = imag(V);
    
    % Jfull = J2n
    Zfull = spalloc(M,M,0);
    Efull = speye(M);
    Jfull = [Zfull Efull; -Efull Zfull]; % skew-symmetric matrix J_2n in Hamiltonian formulation
    
    % Jreduced = J2k
    Zred = spalloc(k,k,0);
    Ered = speye(k);
    Jred = [Zred Ered; -Ered Zred]; % skew-symmetric matrix J_2n in Hamiltonian formulation
    
    A     = [phi -psi; psi phi]; 
    Aplus = Jred'*A'*Jfull; % the symplectic inverse; Aplus*A = I
    
    options.model.A     = A;
    options.model.J     = Jred;
    options.model.Morig = Morig;
    
    figure
    semilogy(diag(Sigma));
    
    u_start = Aplus*u_start;
    
    % dimension of resulting reduced system
    M = 2*k;
    
elseif (strcmp(ROM_basis,'POD'))
    
    Morig = M; % number of nodes
    
    S = snapshots.u_plot(:,1:skip:end);

    % compute the economic SVD of S: S = K*Sigma*L'
    % K is Ndim x Ndim, Sigma is Ndim x Ndim, L is Nt x Ndim
    % in this case Ndim<<Nt
    [K,Sigma,L] = svd(S);

    % construct reduced basis from K, such that V'*V = I
    V = K(:,1:k);

    options.model.V = V;
    options.model.Morig = Morig;
    
    % transformed initial condition
    u_start = V'*u_start;

    % number of unknowns
    M = k;
    
    figure
    semilogy(diag(Sigma));
    
else
    
    error('ROM_basis not available');
end

options.model.ROM_basis = ROM_basis;



