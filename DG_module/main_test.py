import numpy as np
import matplotlib.pyplot as plt
import pdb
import matplotlib.animation as animation
import DG_routines
import DG_solver
import DG_models

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

    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    writergif = animation.PillowWriter(fps=5)

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(sol_list), interval=20, blit=True)

    # save the animation as mp4 video file
    anim.save(gif_name + '.gif',writer=writergif)

K = 20
N = 20

velocity = .3
inflow = 1.
inflow_noise = 0.01

params = {'velocity': velocity,
          'inflow': inflow,
          'inflow_noise': inflow_noise}
error = []
N_vec = [5]
for N in N_vec:
    advection_model = DG_models.Advection(xmin=0,xmax=1,K=K,N=N,
                                   integrator='BDF2',
                                   params=params,
                                   stabilizer_type='filter',Nc=2,s=36
                                   )

    xVec = np.reshape(advection_model.x, (N + 1) * K, 'F')

    FinalTime = 15.

    #mu = 0.4
    #sigma = .075
    #uinit_func = lambda x:1/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*((x-mu)/sigma)**2) + inflow
    uinit_func = lambda x: inflow*np.ones(x.shape)
    uinit = uinit_func(xVec)


    u, time = advection_model.solve(uinit, t_end=FinalTime, step_size=0.1)

    '''
    true_end = uinit_func(xVec - params['velocity'] * time[-1])

    error.append(np.linalg.norm(u[-1] - true_end) / np.linalg.norm(true_end))
    '''
plt.figure()
plt.plot(xVec,uinit,linewidth=2,label='u init')
plt.plot(xVec,u[-1],linewidth=2,label='u approx')
plt.show()

animateSolution(xVec,time[0:-1:2],u[0:-1:2],gif_name='adv_simulation')






'''
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(xVec,uinit.flatten('F'),linewidth=2,label='u init')
ax1.plot(xVec,u[-1],linewidth=2,label='u approx')
ax1.plot(xVec,true_end,'--',linewidth=2,label='u true')
ax1.grid(True)
ax1.legend()

ax2.loglog(N_vec,error,'.-',linewidth=2,markersize=15)
ax2.grid(True)

plt.show()
'''
