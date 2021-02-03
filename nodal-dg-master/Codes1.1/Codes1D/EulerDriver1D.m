% Driver script for solving the 1D Euler equations
Globals1D;

% Polynomial order used for approximation 
N = 1;
K = 250
% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0.0, 1.0, K);

% Initialize solver and construct grid and metric
StartUp1D;
gamma = 1.4;

% Set up initial conditions -- Sod's problem
MassMatrix = inv(V')/V;
cx = ones(Np,1)*sum(MassMatrix*x,1)/2; 

rho = ones(Np,K).*( (cx<0.5) + 0.125*(cx>=0.5));
%rho = exp(-0.5*((x-0.5)/0.1).^2) + 1;
rhou = zeros(Np,K);
Ener = ones(Np,K).*((cx<0.5) + 0.1*(cx>=0.5))/(gamma-1.0);
%Ener = exp(-0.5*((x-0.5)/0.1).^2) + 0.2;
FinalTime = 0.2;

% Solve Problem
[rho,rhou,Ener] = Euler1D(rho,rhou,Ener,FinalTime);

%%
rhoVec = reshape(rho,(N+1)*K,1);
rhouVec = reshape(rhou,(N+1)*K,1);
EnerVec = reshape(Ener,(N+1)*K,1);

plot(rhoVec)
hold on
plot(rhouVec./rhoVec)
legend('rho','u')