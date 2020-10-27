import numpy as np
import matplotlib.pyplot as plt
import pdb
import pipe_flow_global as pipe_flow
import pipe_flow_adjoint as pipe_flow_adjoint
import matplotlib.animation as animation
import scipy.integrate as int
import scipy.special as spec

N = 3
K = 100

xmin = 0.
xmax = 5000.

integrator = True

poly = 'legendre'

c = 1227
rho0 = 870
p0 = 5e5
pamb = 1e5

diameter = 0.508

DG_model_pipe = pipe_flow.Pipe1D(xmin=xmin,xmax=xmax,K=K,N=N,c=c,rho0=rho0,p0=p0, diameter=diameter,poly=poly)
DG_model_pipe.StartUp()

xVec = np.reshape(DG_model_pipe.x, (N + 1) * K, 'F')
xl = 1548
tl = np.array([[20.,1000.]])

Cv = 4.05e-4
mu = 1.04e-1
initInflow = 2.
initOutPres = 5e5

FinalTime = 30.

pinit = initOutPres*np.ones(DG_model_pipe.x.shape)
q1init = ((pinit - DG_model_pipe.p0) / (DG_model_pipe.c ** 2) + DG_model_pipe.rho0) * DG_model_pipe.A
q2init = initInflow*np.ones((N+1,K))*q1init

solq1,solq2, time = DG_model_pipe.solve(q1init,q2init, FinalTime=FinalTime,implicit=integrator,stepsize=2e0,
                                        xl=xl,tl=tl,leak_type='discharge',Cv=Cv,pamb=pamb,mu=mu,
                                        initInflow=initInflow,initOutPres=initOutPres)
#%%
rhoA = []
u = []
p = []
pBar = []
rho = []
for i in range(len(solq1)):
    rhoA.append(solq1[i].flatten('F'))
    rho.append((solq1[i]/DG_model_pipe.A).flatten('F'))
    u.append(np.divide(solq2[i],solq1[i]).flatten('F'))
    p.append(np.reshape(c*c*(solq1[i]/DG_model_pipe.A-rho0)+p0, (N + 1) * K, 'F'))
    pBar.append(1e-5*np.reshape(c*c*(solq1[i]/DG_model_pipe.A-rho0)+p0, (N + 1) * K, 'F'))


DG_model_pipe_adjoint = pipe_flow_adjoint.Pipe1D_Adjoint(DG_model_pipe)

obs = np.zeros(2*(N+1)*K)

lol = DG_model_pipe_adjoint.SolveAdjoint(solq1[-1],solq2[-1],obs)