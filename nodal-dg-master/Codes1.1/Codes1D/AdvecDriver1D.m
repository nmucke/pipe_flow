% Driver script for solving the 1D advection equations
Globals1D;

% Order of polymomials used for approximation 
N = 4;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0.0,2*pi,2);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
uInit = sin(x);

% Solve Problem
FinalTime = .1;
[u] = Advec1D(uInit,FinalTime);

truesol = sin(x-2*pi*FinalTime);

plot(x(:),u(:),'linewidth',2.0)
hold on
plot(x(:),truesol(:),'--','linewidth',2.0)
plot(x(:),uInit(:),'--','linewidth',2.0)
grid on
