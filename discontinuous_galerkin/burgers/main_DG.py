import numpy as np
import matplotlib.pyplot as plt
import pdb
import DG_routines as DG
import matplotlib.animation as animation
import scipy.integrate as int
import scipy.special as spec
from scipy.optimize import newton

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
poly = 'legendre'
integrator = False
error = []
N_vec = [2]
for N in N_vec:
    #N = 2
    K = 100

    xmin = 0
    xmax = 1.

    burgers = DG.Burgers1D(xmin=xmin,xmax=xmax,K=K,N=N,poly=poly)
    burgers.StartUp()

    xVec = np.reshape(burgers.x, (N + 1) * K, 'F')

    FinalTime = 1.8

    ts = 2.
    eps = 0.0000000001*np.pi
    #uinit = -np.tanh((burgers.x+0.5)/eps) + 1
    uinit = 1/(2*np.pi*ts)*np.sin(2*np.pi*burgers.x)
    #uinit = np.sin(burgers.x)

    u, time = burgers.solve(uinit, FinalTime=FinalTime,implicit=False,stepsize=2e0)

    u_mat = np.asarray(u)

    true_end_func = lambda y: 1/(2*np.pi*ts)*np.sin(2*np.pi*(burgers.x.flatten('F')-y*time[-1]))-y
    true_end = newton(true_end_func,uinit.flatten('F'),tol=1e-12)
    #true_end=-np.tanh((burgers.x.flatten('F')+0.5-time[-1])/eps) + 1
    #true_end = np.sin(burgers.x.flatten('F')-time[-1])

    error.append(np.linalg.norm(u[-1]-true_end)/np.linalg.norm(true_end))

plt.figure()
plt.plot(xVec,uinit.flatten('F'),linewidth=2,label='u init')
plt.plot(xVec,u[np.int(len(u)/2)],linewidth=2,label='u middle')
plt.plot(xVec,u[-1],linewidth=2,label='u end')
plt.grid(True)
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('u')
plt.show()

'''
plt.figure()
plt.plot(xVec,uinit.flatten('F'),linewidth=2,label='u init')
plt.plot(xVec,u[-1],linewidth=2,label='u approx')
plt.plot(xVec,true_end,'--',linewidth=2,label='u true')
plt.grid(True)
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('u')
plt.show()

plt.figure()
plt.loglog(N_vec,error,'.-',linewidth=2,markersize=15)
plt.grid(True)
plt.legend()
plt.xlabel('N')
plt.ylabel('error')
plt.show()
'''

#animateSolution(xVec,time[0:-1],u[0:-1])


