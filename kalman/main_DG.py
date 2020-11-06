import numpy as np
import matplotlib.pyplot as plt
import pdb
import pipe_flow as pipe_flow
import pipe_flow_kalman as pipe_flow_kalman

import matplotlib.animation as animation
import scipy.integrate as int
import scipy.special as spec

# animation function.  This is called sequentially
def animateSolution(x,time,sol_list,gif_name='pipe_flow_simulation'):
    fig = plt.figure()
    ax = plt.axes(xlim=(x[0], x[-1]), ylim=(np.min(sol_list),np.max(sol_list)))
    #ax = plt.axes(xlim=(x[0], x[-1]), ylim=(-1,1))

    #ax = plt.axes(xlim=(x[0], x[-1]), ylim=(np.min(sol_list),1003))

    plt.grid(True)
    line, = ax.plot([], [], lw=2)

    # initialization function: plot the background of each frame
    def init():
        line.set_data([], [])
        return line,

    def animate(i):
        plt.title(str(time[i]))
        y = sol_list[i]
        line.set_data(x, y)
        return line,

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(sol_list), interval=20, blit=True)

    # save the animation as mp4 video file
    anim.save(gif_name + '.mp4',writer=writer)

poly = 'legendre'
implicit = False
N = 1
K = 100

xmin = 0.
xmax = 5000.
c = 1227#1400
rho0 = 870#1000
p0 = 5e5
pamb = 1e5
diameter = 0.508
mu = 1.04e-1
initInflow = 2.
initOutPres = 5e5

FinalTime = 100

observation_model = pipe_flow.Pipe1D(xmin=xmin,xmax=xmax,K=K,N=N,c=c,rho0=rho0,p0=p0,
                          diameter=diameter,poly=poly)
observation_model.StartUp()
xVec = np.reshape(observation_model.x, (N + 1) * K, 'F')

xl_true = 4548
Cv_true = 4.05e-4
tl = np.array([[20.,1000.]])

pinit = p0*np.ones(observation_model.x.shape)
q1init = ((pinit - observation_model.p0) / (observation_model.c ** 2) + observation_model.rho0) * observation_model.A
q2init = initInflow*np.ones((N+1,K))*q1init

solq1,solq2, time = observation_model.solve(q1init,q2init, FinalTime=FinalTime,implicit=implicit,stepsize=2e0,
                                        xl=xl_true,tl=tl,leak_type='discharge',Cv=Cv_true,pamb=pamb,mu=mu,initInflow=initInflow,initOutPres=initOutPres)
#%%
animateSolution(xVec,time,solq2)
def observation_operator(q1,q2):
    u = []
    for i in range(len(solq1)):
        u.append(np.divide(q2[i],q1[i]).flatten('F'))
    u = np.asarray(u)
    obs = u[:,-1]
    return obs

observations = observation_operator(solq1,solq2)

plt.figure()
plt.plot(time,observations)
plt.grid()
plt.show()

obs_noise = 2
R = np.array([obs_noise])

N_ensembles = 20


observation_model = pipe_flow_kalman.Pipe1D(xmin=xmin,xmax=xmax,K=K,N=N,c=c,rho0=rho0,p0=p0,
                          diameter=diameter,poly=poly)
observation_model.StartUp()
xVec = np.reshape(observation_model.x, (N + 1) * K, 'F')

xl_true = 4548
Cv_true = 4.05e-4
tl = np.array([[20.,1000.]])

pinit = p0*np.ones(observation_model.x.shape)
q1init = ((pinit - observation_model.p0) / (observation_model.c ** 2) + observation_model.rho0) * observation_model.A
q2init = initInflow*np.ones((N+1,K))*q1init

solq1,solq2, time = observation_model.solve(q1init,q2init, FinalTime=FinalTime,implicit=implicit,stepsize=2e0,
                                        xl=xl_true,tl=tl,leak_type='discharge',Cv=Cv_true,pamb=pamb,mu=mu,initInflow=initInflow,initOutPres=initOutPres)
#%%





























