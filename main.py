import numpy as np
import matplotlib.pyplot as plt
import pdb
import DG_routines as DG
import scipy as sci
import scipy.sparse as sps
import matplotlib.animation as animation


import scipy.io


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


alpha = 0
beta = 0

N = 2
#K = 50

xmin = 0.
xmax = 100.

error = []
test_vec = [1000]


for K in test_vec:

    c = 1400
    rho0 = 1000
    p0 = 1e5

    DG_model_pipe = DG.DG_1D(xmin=xmin, xmax=xmax, K=K, N=N)
    DG_model_pipe.StartUp()

    xVec = np.reshape(DG_model_pipe.x, (N + 1) * K, 'F')

    xl = 50

    mu1 = 50.
    mu2 = 60
    sigma = .01
    #q1init =  3 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * np.power((DG_model_pipe.x - mu1) / sigma, 2))
    #q1init +=  3 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * np.power((DG_model_pipe.x - mu2) / sigma, 2))
    #q1init += np.random.normal(0,1,q1init.shape)
    #q1init += rho0


    #q2init = np.zeros((N+1,K))

    #q1init = np.zeros(DG_model_pipe.x.shape)
    #q1init[np.argwhere(DG_model_pipe.x>=20)[:,0],np.argwhere(DG_model_pipe.x>=20)[:,1]] =rho0 + 1.25
    #q1init[np.argwhere(DG_model_pipe.x<20)[:,0],np.argwhere(DG_model_pipe.x<20)[:,1]] = rho0 + 3.5
    q2init = np.zeros((N+1,K))
    q1init = rho0*np.ones((N+1,K))

    solq1,solq2, time = DG_model_pipe.Pipe1D(q1init,q2init, FinalTime=0.2,c=c,rho0=rho0,p0=p0,xl=xl)
#%%
rho = []
u = []
p = []
for i in range(len(solq1)):
    rho.append(np.reshape(solq1[i], (N + 1) * K, 'F'))
    u.append(np.reshape(np.divide(solq2[i],solq1[i]), (N + 1) * K, 'F'))
    p.append(np.reshape(c*c*(solq1[i]-rho0)+p0, (N + 1) * K, 'F'))

#%%

xVec = np.reshape(DG_model_pipe.x, (N + 1) * K, 'F')
plt.figure()
plt.plot(xVec,u[-1])
plt.grid(True)
plt.legend(['u'])
plt.show()

plt.figure()
plt.plot(xVec,rho[-1])
plt.grid(True)
plt.legend(['rho'])
plt.show()

plt.figure()
plt.plot(xVec,p[-1])
plt.grid(True)
plt.legend(['p'])
plt.show()


#%%
#animateSolution(xVec,time[0:-1:20],u[0:-1:20])

