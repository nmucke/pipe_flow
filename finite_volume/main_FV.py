import numpy as np
import matplotlib.pyplot as plt
import pdb
import FV_routines as FV
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
    anim.save(gif_name + '.mp4',writer=writer)

'''
c = 1400
rho0 = 1000
p0 = 1e5
diameter = 0.5

xmin, xmax = 0,2*np.pi
num_volumes = 1000

error = []
vec = [10,50,100,150,200,300,400,500,1000,2000,3000,4000,5000,6000]
for num_volumes in vec:


    advection_model = FV.Advection1D(xmin=xmin,xmax=xmax,num_volumes=num_volumes,c=1)


    FinalTime = 10
    u0 = np.sin(advection_model.xmid)
    sol = advection_model.solve(u0,FinalTime=FinalTime)

    true_sol = lambda t: np.sin(advection_model.xmid-1*t)
    pdb.set_trace()

    error.append(1/np.sqrt(num_volumes)*np.linalg.norm(true_sol(FinalTime)-sol[-1]))


plt.figure()
plt.plot(sol[-1])
plt.plot(true_sol(10.2))
plt.show()

plt.figure()
plt.loglog(vec,error)
plt.loglog(vec,np.divide(1,vec))
plt.grid()
plt.show()
'''


xmin, xmax = 0,100
num_volumes = 1000

error = []
vec = [1000]
for num_volumes in vec:


    LSWE_model = FV.LSWE(xmin=xmin,xmax=xmax,num_volumes=num_volumes,d0=1,g=9.81)

    FinalTime = 0.25

    mu = xmax/2
    sigma = 0.1
    q1init = 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * np.power((LSWE_model.xmid - mu) / sigma, 2))
    #q1init = np.cos(LSWE_model.xmid)

    q2init = np.zeros(LSWE_model.x.shape)
    q2init = q2init[0:-1]
    solq1,solq2,time = LSWE_model.solve_RK(q1init,q2init,FinalTime=FinalTime)

    #true_sol = lambda t: np.sin(LSWE_model.xmid-1*t)

    #error.append(1/np.sqrt(num_volumes)*np.linalg.norm(true_sol(FinalTime)-sol[-1]))


plt.figure()
plt.plot(solq1[-1])
plt.plot(solq2[-1])
#plt.plot(true_sol(10.2))
plt.show()

animateSolution(LSWE_model.xmid,time[0:-1:20],solq2[0:-1:20],gif_name='pipe_flow_simulation')
'''
plt.figure()
plt.loglog(vec,error)
plt.loglog(vec,np.divide(1,vec))
plt.grid()
plt.show()
'''