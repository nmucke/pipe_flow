import numpy as np
import matplotlib.pyplot as plt
import pdb
import pipe_flow_global as pipe_flow
import pipe_flow_adjoint as pipe_flow_adjoint
import matplotlib.animation as animation
import scipy.integrate as int
import scipy.special as spec
import scipy.optimize as opt

N = 1
K = 300

xmin = 0.
xmax = 5000.

implicit = False

poly = 'legendre'

c = 1227
rho0 = 870
p0 = 5e5
pamb = 1e5

diameter = 0.508

Cv = 4.05e-4
mu = 1.04e-1
initInflow = 2.
initOutPres = 5e5

FinalTime =5.
tl = np.array([[0.,1000.]])

pipe_obs = pipe_flow.Pipe1D(xmin=xmin,xmax=xmax,K=K,N=N,c=c,rho0=rho0,p0=p0, diameter=diameter,poly=poly)
pipe_obs.StartUp()

xVec = np.reshape(pipe_obs.x, (N + 1) * K, 'F')
xl_true = 4548

pinit = initOutPres*np.ones(pipe_obs.x.shape)
q1init = ((pinit - pipe_obs.p0) / (pipe_obs.c ** 2) + pipe_obs.rho0) * pipe_obs.A
q2init = initInflow*np.ones((N+1,K))*q1init

solq1_obs,solq2_obs, time = pipe_obs.solve(q1init,q2init, FinalTime=FinalTime,implicit=implicit,stepsize=2e0,
                                        xl=xl_true,tl=tl,leak_type='discharge',Cv=Cv,pamb=pamb,mu=mu,
                                        initInflow=initInflow,initOutPres=initOutPres)


obs = np.zeros(solq1_obs.shape)
obs[:,-1] = solq2_obs[:,-1]/solq1_obs[:,-1]
obs = np.concatenate((obs,obs),axis=1)

pipe_solve = pipe_flow.Pipe1D(xmin=xmin,xmax=xmax,K=K,N=N,c=c,rho0=rho0,p0=p0, diameter=diameter,poly=poly)
pipe_solve.StartUp()

pinit = initOutPres*np.ones(pipe_solve.x.shape)
q1init = ((pinit - pipe_solve.p0) / (pipe_solve.c ** 2) + pipe_solve.rho0) * pipe_obs.A
q2init = initInflow*np.ones((N+1,K))*q1init

xVec = np.reshape(pipe_solve.x, (N + 1) * K, 'F')
xl_0 = 2500
xl = xl_0
error = 1e3
reg = 1e-6

def J(q1,q2,obs,xl,xl0,reg):

    cost = 1/len(q1) * np.sum(np.power(q2[:,-1]/q1[:,-1]-obs[:,-1],2)) + reg/2*np.abs(xl-xl0)**2

    return cost

def cost_newton(xl):

    print(xl)

    xl = pipe_solve.deltax + xl

    solq1, solq2, time = pipe_solve.solve(q1init, q2init, FinalTime=FinalTime, implicit=implicit, stepsize=2e0,
                                          xl=xl, tl=tl, leak_type='discharge', Cv=Cv, pamb=pamb, mu=mu,
                                          initInflow=initInflow, initOutPres=initOutPres)

    cost = np.sum(np.power(solq2[:,[3,10,40,56,100,137,250,-1]]/solq1[:,[3,10,40,56,100,137,250,-1]]-obs[:,[3,10,40,56,100,137,250,-1]],2)) + reg/2*np.abs(xl-2500)**2

    print(cost)
    return cost

xxx = opt.minimize(cost_newton,np.array([2500.]),tol=1e-8)
pdb.set_trace()
while error > 1e-6:

    solq1,solq2, time = pipe_solve.solve(q1init,q2init, FinalTime=FinalTime,implicit=implicit,stepsize=2e0,
                                            xl=xl,tl=tl,leak_type='discharge',Cv=Cv,pamb=pamb,mu=mu,
                                            initInflow=initInflow,initOutPres=initOutPres)

    pres = c*c*(solq1/pipe_solve.A-rho0)+p0
    rho = solq1/pipe_solve.A
    leakage = []
    for i in range(len(time)):
        l = pipe_solve.Leakage(time[i],pipe_solve.xElementL,tl,pressure=pres[i],rho=rho[i]).flatten('F')
        leakage.append(pipe_solve.invM_global.dot(l))
    leakage = np.asarray(leakage)
    leakage = np.concatenate((leakage,np.zeros(leakage.shape)),axis=1)

    adjoint_sol = pipe_solve.SolveAdjoint(solq1, solq1, obs, time)
    adj1 = adjoint_sol[:,0:(N+1)*K]
    adj2 = adjoint_sol[:,-((N+1)*K):]

    adjDotLeakage = []
    for i in range(len(time)):
        adjDotLeakage.append(np.dot(adjoint_sol[-i],1e8*leakage[i]))

    xl_step = (reg*np.abs(xl_0-xl) + 2*np.sum(adjDotLeakage))
    pdb.set_trace()
    xl += xl_step

    error = np.abs(xl_step)

plt.figure()
plt.plot(xVec,adj2[1,:]/adj1[1,:])
plt.plot(xVec,adj2[50,:]/adj1[50,:])
plt.plot(xVec,adj2[-1,:]/adj1[-1,:])
plt.grid()
plt.show()





#adjoint_sol = DG_model_pipe.SolveAdjoint(solq1, solq1, obs,time)
'''
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

leaked_mass = []
for t in range(len(time)):
    leak = DG_model_pipe.Leakage(time[t],DG_model_pipe.xElementL,tl,pressure=p[t],rho=rho[t])
    leaked_mass.append(leak[0,DG_model_pipe.xElementL])

pressure_leak = []
for t in range(len(time)):
    mean_pressure = np.mean(np.reshape(p[t], (DG_model_pipe.Np, DG_model_pipe.K), 'F')[:, DG_model_pipe.xElementL])
    pressure_leak.append(mean_pressure)
benjamin_data = np.genfromtxt('benjamin_data.csv',delimiter=',')
benjamin_time = benjamin_data[:,0]
benjamin_mass_leak = benjamin_data[:,1]

plt.figure()
plt.plot(time,leaked_mass,label='DG',linewidth=2)
plt.plot(benjamin_time,benjamin_mass_leak,label='Rosa',linewidth=2)
plt.grid(True)
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Mass Flow Through Leak (kg/s)')
plt.savefig('leaked_mass')
plt.show()

#DG_model_pipe_adjoint = pipe_flow_adjoint.Pipe1D_Adjoint(DG_model_pipe)

#obs = np.zeros(2*(N+1)*K)

#lol = DG_model_pipe_adjoint.SolveAdjoint(solq1[-1],solq2[-1],obs)

'''