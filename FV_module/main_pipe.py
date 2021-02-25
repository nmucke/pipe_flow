import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy.io import loadmat
import matplotlib.animation as animation

import FV_models

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

    writergif = animation.PillowWriter(fps=5)

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(sol_list), interval=20, blit=True)

    # save the animation as mp4 video file
    anim.save(gif_name + '.gif', writer=writergif)

params = {'diameter': 0.508,
          'rho0': 870,
          'pamb': 1e5,
          'c': 1227,
          'p0': 5e5}

num_volumes = 300


xmin = 0
xmax = 50
Pipe = FV_models.Pipe(xmin=xmin,xmax=xmax,num_volumes=num_volumes,
                      params=params,integrator='LowStorageRK')


A = params['diameter']**2/4 * np.pi
mu = xmax/2
sigma = 1
q1_init_func = lambda x: (1 + 1/(sigma*np.sqrt(2*np.pi))
                         *np.exp(-0.5*((x-mu)/sigma)**2)) * A
q1_init = q1_init_func(Pipe.xmid)

q2_init = np.zeros(Pipe.x.shape)

sol, time = Pipe.solve(q1_init,q2_init,t_end=0.1,step_size=0.001)

u = []
rho = []

for i in range(len(time)):
    rho.append(sol[i][0:num_volumes]/A)
    u.append(sol[i][-(num_volumes+1):]/np.dot(Pipe.I_p,rho[i])/A)


animateSolution(Pipe.x, time[0:-1:5], u[0:-1:5])

plt.figure()
plt.plot(Pipe.xmid,rho[-1],'-',linewidth=1.5)
#plt.plot(Pipe.x,u,'-',linewidth=1.5)
plt.show()