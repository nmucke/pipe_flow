close all
clear all 
clc 

Globals1D;

% Order of polymomials used for approximation 
N = 2; 

% Gravitational acceleration
g = 9.81;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0,100,500); 

% Initialize solver and construct grid and metric 
StartUp1D; 

global p0 rho0 c
c = 1400;
rho0 = 1000;
p0 = 1e5;

% Set initial conditions 
mu = 50;
sigma = 1;
h =  1 / (sigma * sqrt(2 * pi)) * exp(-0.5 * ((x - mu) / sigma).^2);
h = h + rho0;
hu = zeros(size(x)); 

dt = 1e0;

% Solve Problem 
FinalTime = 10; 
[h,hu] = pipe_flow_implicit(h,hu,FinalTime,dt); 





