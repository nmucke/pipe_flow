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
for K in [100,200,300,400]:
    N = 3
    #K = 500

    xmin = 0.
    xmax = 100.

    c = 1400
    rho0 = 1000
    p0 = 1e5

    diameter = 0.03568248232

    DG_model_pipe = DG.Pipe1D(xmin=xmin,xmax=xmax,K=K,N=N,c=c,rho0=rho0,p0=p0, diameter=diameter)
    DG_model_pipe.StartUp()

    xVec = np.reshape(DG_model_pipe.x, (N + 1) * K, 'F')
    xl = 50
    tl = np.array([[.4,.5]])

    mu1 = 80.
    sigma = .01
    q1init =  1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * np.power((DG_model_pipe.x - mu1) / sigma, 2)) + rho0
    q1init *= DG_model_pipe.A
    #q1init = rho0*np.ones((N+1,K))*DG_model_pipe.A

    q2init = np.zeros((N+1,K))

    solq1,solq2, time = DG_model_pipe.solve(q1init,q2init, FinalTime=.2,xl=xl,tl=tl,implicit=False,stepsize=1e-3)
    #%%
    rhoA = []
    u = []
    p = []
    for i in range(len(solq1)):
        rhoA.append(np.reshape(solq1[i], (N + 1) * K, 'F'))
        u.append(np.reshape(np.divide(solq2[i],solq1[i]), (N + 1) * K, 'F'))
        p.append(np.reshape(c*c*(solq1[i]-rho0)+p0, (N + 1) * K, 'F'))

    '''
    initial_total_mass = 0
    end_total_mass = 0
    for i in range(K):
        initial_total_mass += int.simps(solq1[0][:,i],DG_model_pipe.x[:,i])
        end_total_mass += int.simps(solq1[-1][:,i],DG_model_pipe.x[:,i])


    print('Initial total mass: ' + str(initial_total_mass))
    print('End total mass: ' + str(end_total_mass))
    print('Mass difference: ' + str(initial_total_mass-end_total_mass))
    '''
    xVec = np.reshape(DG_model_pipe.x, (N + 1) * K, 'F')

    plt.plot(xVec, u[-1],label=str(N))
    plt.grid(True)
plt.legend()
plt.show()
    #%%
'''
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
'''

#%%
animateSolution(xVec,time[0:-1:20],u[0:-1:20])

