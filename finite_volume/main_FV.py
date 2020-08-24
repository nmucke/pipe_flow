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


xmin=0
xmax=10
num_volumes=100

d0 = 1
g = 9.81
LSWE_model = FV.LSWE(xmin=xmin,xmax=xmax,num_volumes=num_volumes,d0=d0,g=g)

mu = 5
sigma = 0.05

q1init = 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-0.5*np.power((LSWE_model.xmid-mu)/sigma,2))
q2init = np.zeros(LSWE_model.x.shape)

solq1,solq2, time = LSWE_model.solve_RK(q1init,q2init, FinalTime=0.5)
solq1_array = np.asarray(solq1)
solq2_array = np.asarray(solq2)

plt.figure()
plt.plot(solq2_array[-1])
plt.show()

skip_idx = np.int(np.round(len(solq1)/100))
animateSolution(LSWE_model.x,time[0:-1:skip_idx],solq2[0:-1:skip_idx],gif_name='LSWE_simulation')




















'''
xmin=-1
xmax=1
num_volumes=1000
c=1
test_vec = [10]
#test_vec = [100]
error = []
for num_volumes in test_vec:
    adv_model = FV.Burgers1D(xmin=xmin,xmax=xmax,num_volumes=num_volumes,c=c)
    u0 = 1/(0.5*np.sqrt(2*np.pi)) * np.exp(-0.5*np.power((adv_model.xmid-np.pi)/0.5,2))
    u0 = adv_model.xmid
    sol, time = adv_model.solve_RK(u0, FinalTime=10)
    sol_array = np.asarray(sol)


    true_sol = lambda t: adv_model.xmid/(t+1)
    true_sol_array = []
    for i in range(len(time)):
        true_sol_array.append(true_sol(time[i]))
    true_sol_array = np.asarray(true_sol_array)

    e = np.mean(1/np.sqrt(num_volumes)*np.linalg.norm(sol_array-true_sol_array,axis=0))
    error.append(e)



plt.figure()
plt.loglog(test_vec,error,'.-',linewidth=2,markersize=15,label='FV Error')
plt.loglog(test_vec,np.divide(1,test_vec),'.-',linewidth=2,markersize=15,label='h^-1')
plt.grid()
plt.legend()
plt.show()


plt.figure()
plt.plot(sol[-1])
plt.plot(true_sol_array[-1])
plt.grid()
plt.show()

skip_idx = np.int(np.round(len(sol)/100))
animateSolution(adv_model.xmid,time[0:-1:skip_idx],sol[0:-1:skip_idx],gif_name='advection_simulation')

'''