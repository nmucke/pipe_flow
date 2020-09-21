import numpy as np
import matplotlib.pyplot as plt
import pdb
import DG_routines as DG
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
plt.figure()
solsu = []
solsrhoa = []
solsp = []
for N in [2]:
    #N = 2
    K = 250

    xmin = 0.
    xmax = 5000.

    c = 1227#1400
    rho0 = 870#1000
    p0 = 5e5
    pamb = 1e5

    diameter = 0.508

    DG_model_pipe = DG.Pipe1D(xmin=xmin,xmax=xmax,K=K,N=N,c=c,rho0=rho0,p0=p0, diameter=diameter)
    DG_model_pipe.StartUp()

    xVec = np.reshape(DG_model_pipe.x, (N + 1) * K, 'F')
    xl = 1000
    tl = np.array([[20.,1000.]])
    Cv = 1.04e-4
    mu = 1.04e-1
    initInflow = 2.
    initOutPres = 5e5

    FinalTime = 20.

    mu1 = xmax/2
    sigma = .5
    #q1init =  100 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * np.power((DG_model_pipe.x - mu1) / sigma, 2)) + rho0
    #q1init *= DG_model_pipe.A
    pinit = pamb*np.ones(DG_model_pipe.x.shape)
    q1init = ((pinit - DG_model_pipe.p0) / (DG_model_pipe.c ** 2) + DG_model_pipe.rho0) * DG_model_pipe.A#rho0*np.ones((N+1,K))*DG_model_pipe.A
    q2init = initInflow*np.ones((N+1,K))

    solq1,solq2, time = DG_model_pipe.solve(q1init,q2init, FinalTime=FinalTime,implicit=True,stepsize=1e-1,
                                            xl=xl,tl=tl,leak_type='discharge',Cv=Cv,pamb=pamb,mu=mu,initInflow=initInflow,initOutPres=initOutPres)
    #%%
    rhoA = []
    u = []
    p = []
    rho = []
    for i in range(len(solq1)):
        rhoA.append(solq1[i].flatten('F'))
        rho.append((solq1[i]/DG_model_pipe.A).flatten('F'))
        u.append(np.divide(solq2[i],solq1[i]).flatten('F'))
        p.append(np.reshape(c*c*(solq1[i]/DG_model_pipe.A-rho0)+p0, (N + 1) * K, 'F'))

    solsu.append(u)
    solsrhoa.append(rhoA)
    solsp.append(p)

    initial_total_mass = 0
    end_total_mass = 0
    for i in range(K):
        initial_total_mass += int.simps(np.reshape(solq1[0],(N+1,K),'F')[:,i],DG_model_pipe.x[:,i])
        end_total_mass += int.simps(np.reshape(solq1[-1],(N+1,K),'F')[:,i],DG_model_pipe.x[:,i])

    total_mass = []
    for t in range(len(time)):
        mass = 0
        for i in range(K):
            mass += int.simps(np.reshape(solq1[t],(N+1,K),'F')[:,i],DG_model_pipe.x[:,i])
        total_mass.append(mass)


    print('Initial total mass: ' + str(initial_total_mass))
    print('End total mass: ' + str(end_total_mass))
    print('Mass difference: ' + str(initial_total_mass-end_total_mass))

    xVec = DG_model_pipe.x.flatten('F')
    plt.plot(xVec, u[-1],label=str(K))
plt.grid(True)
plt.legend()
plt.show()

plt.figure()
plt.plot(time,total_mass)
plt.grid(True)
plt.legend(['Total Mass'])
plt.show()
    #%%

xVec = np.reshape(DG_model_pipe.x, (N + 1) * K, 'F')
plt.figure()
plt.plot(xVec,u[-1])
plt.grid(True)
plt.legend(['u'])
plt.show()

plt.figure()
plt.plot(xVec,rhoA[-1])
plt.grid(True)
plt.legend(['rhoA'])
plt.show()

plt.figure()
plt.plot(xVec,p[-1])
plt.grid(True)
plt.legend(['p'])
plt.show()


#%%
animateSolution(xVec,time[0:-1],u[0:-1])

