close all
clear all 
clc 
% Driver script for solving the 1D non?linear SWE equations 
%
% By Allan P. Engsig-Karup, apek@imm.dtu.dk.
Globals1D;

% Order of polymomials used for approximation 
N = 2; 

% Gravitational acceleration
g = 9.81;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0,50,201); 

% Initialize solver and construct grid and metric 
StartUp1D; 

% Inital water levels 
hl = 3.5; 
hr = 1.25; 

% Set initial conditions 
h = zeros(size(x)); 
h(x<20) = hl; 
h(x>20) = hr; 
hu = zeros(size(x)); 

%mu = 0.;
%sigma=0.1;
%h = 1+ 10 / (sigma * sqrt(2 * pi)) * exp(-0.5 * ((x - mu) / sigma).^2);

dt = 0.01;

% Solve Problem 
FinalTime = 2.5; 
[h,hu] = NSWE1D_implicit(h,hu,FinalTime,g,dt); 

% Load data for reference solution
load('AssignmentIData/dambreakdata.mat')

subplot(2,1,1)
hold on
plot(x,h,'r')
%ylim([1 4])
legend('Computed','Exact')

subplot(2,1,2)
hold on 
plot(x,hu,'r')
%ylim([-5 10])
legend('Computed','Exact')
%print ?f1 ?dpng dambreakN1.png





