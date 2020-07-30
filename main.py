import numpy as np
import matplotlib.pyplot as plt
import pdb
import DG_routines as DG
import matplotlib.animation as animation

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
    anim.save(gif_name + '.avi',writer=writer)

N = 2
K = 200

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
tl = np.array([[0.05,0.1]])

mu1 = 50.
mu2 = 60
sigma = .01
q1init =  3 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * np.power((DG_model_pipe.x - mu1) / sigma, 2))
q1init += rho0

#q1init = rho0*np.ones((N+1,K))

q2init = np.zeros((N+1,K))

solq1,solq2, time = DG_model_pipe.solve(q1init,q2init, FinalTime=.1,xl=xl,tl=tl,implicit=True,stepsize=1e-2)
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
animateSolution(xVec,time[0:-1:5],u[0:-1:5])

