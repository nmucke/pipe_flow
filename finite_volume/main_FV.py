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
N = 200
A = -2*np.eye(N+1) + np.eye(N+1,k=-1) + np.eye(N+1,k=1)
A[0,0],A[0,1] = 1,0
A[-1,-1],A[-1,-2] = 1,0


x = np.linspace(0,1,N+1)
dx = x[1]-x[0]
initCondition = 6*np.sin(np.pi*x)
k = 1.

f = lambda t:-k*np.pi**2*6*np.sin(np.pi*x)*np.exp(-k*np.pi**2*t)+6*np.pi**2*np.sin(np.pi*x)*np.exp(-k*np.pi**2*t)+36*np.sin(np.pi*x)*np.sin(np.pi*x)*np.exp(-2*k*np.pi**2*t)

def rhs(t,y):

    #y_diff = np.array([y[1],-y[0]])
    #y_diff = np.array([y[1],-9.81/1*np.sin(y[0])])

    y_diff = k/dx/dx*np.dot(A,y) - y*y + f(t)
    y_diff[0], y_diff[-1] = 0,0
    #y_diff = np.concatenate((np.array([0]),y_diff,np.array([0])))
    return y_diff

truesol = lambda t: 6*np.sin(np.pi*x)*np.exp(-k*(np.pi**2)*t)

vec = [1e0,1e-1,1e-2,1e-3]
error = []
#truesol = lambda t: np.array([2*np.cos(t)+3*np.sin(t),-2*np.sin(t) + 3*np.cos(t) ])
#truesol = lambda t: 6*np.sin(np.pi*x)*np.exp(-k*(np.pi**2)*t)
#system = FV.BDF2(rhs, initCondition, t0=0, te=1, stepsize=1e-6)
#t_vec, true = system.solve()
for stepsize in vec:
    system = FV.BDF2(rhs, initCondition, t0=0, te=1, stepsize=stepsize)
    t_vec, solution = system.solve()

    true = truesol(t_vec[-1])
    e = np.linalg.norm(solution[-1]-true)/np.linalg.norm(true)

    error.append(e)

plt.figure()
plt.plot(solution[-1])
plt.plot(true)
plt.show()
#animateSolution(x,t_vec,solution,gif_name='heat')
plt.figure()
plt.loglog(vec,error)
#plt.loglog(vec,np.divide(1,vec))
plt.loglog(vec,10*np.power(vec,2))
plt.grid()
plt.legend(['error','h2'])
plt.show()
'''
'''
c = 1400
rho0 = 1000
p0 = 1e5
diameter = 0.5

xmin, xmax = 0,2*np.pi
num_volumes = 1000

error = []
vec = [200,300,400]
vec = [0.1,0.095,0.09]
for stepsize in vec:
#for num_volumes in vec:


    advection_model = FV.Advection1D(xmin=xmin,xmax=xmax,num_volumes=num_volumes,c=1)


    FinalTime = 1.
    u0 = np.sin(advection_model.xmid)
    sol,time = advection_model.solve_implicit(u0,FinalTime=FinalTime,stepsize=stepsize)
    #sol,time = advection_model.solve_RK(u0,FinalTime=FinalTime)

    true_sol = lambda t: np.sin(advection_model.xmid-1*t)

    error.append(1/np.sqrt(num_volumes)*np.linalg.norm(true_sol(time[-1])-sol[-1]))


plt.figure()
plt.plot(sol[-1])
plt.plot(true_sol(FinalTime))
plt.show()

plt.figure()
plt.loglog(vec,error)
#plt.loglog(vec,np.divide(1,vec))
plt.loglog(vec,vec)
plt.grid()
plt.legend(['error','h1'])
plt.show()
'''

'''
xmin, xmax = 0,5
num_volumes = 1000

error = []
vec = [1000]
for num_volumes in vec:


    LSWE_model = FV.LSWE(xmin=xmin,xmax=xmax,num_volumes=num_volumes,d0=1,g=9.81)

    FinalTime = 1

    mu = 2.5
    sigma = 0.1
    q1init = 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * np.power((LSWE_model.xmid - mu) / sigma, 2))
    #q1init = np.cos(LSWE_model.xmid)

    q2init = np.zeros(LSWE_model.x.shape)
    solq1,solq2,time = LSWE_model.solve_RK(q1init,q2init,FinalTime=FinalTime)

    #true_sol = lambda t: np.sin(LSWE_model.xmid-1*t)

    #error.append(1/np.sqrt(num_volumes)*np.linalg.norm(true_sol(FinalTime)-sol[-1]))

solq1_right_boundary = []
solq1_left_boundary = []
solq2_right_boundary = []
solq2_left_boundary = []
for i in range(len(solq1)):
    solq1_right_boundary.append(solq1[i][-1])
    solq1_left_boundary.append(solq1[i][0])
    solq2_right_boundary.append(solq2[i][-1])
    solq2_left_boundary.append(solq2[i][0])

plt.figure()
plt.plot(solq2_left_boundary)
plt.plot(solq2_right_boundary)
#plt.plot(true_sol(10.2))
plt.show()

plt.figure()
plt.plot(solq1[-1])
plt.plot(solq2[-1])
#plt.plot(true_sol(10.2))
plt.show()

animateSolution(LSWE_model.xmid,time[0:-1:20],solq2[0:-1:20],gif_name='pipe_flow_simulation')

plt.figure()
plt.loglog(vec,error)
plt.loglog(vec,np.divide(1,vec))
plt.grid()
plt.show()
'''


xmin, xmax = 0,10
num_volumes = 2000

c=1400
rho0=1000
p0=1e5
diameter=0.5

pipe_model = FV.Pipe1D(xmin=xmin,xmax=xmax,num_volumes=num_volumes,c=c,rho0=rho0,p0=p0, diameter=diameter)

FinalTime = 0.01

mu = xmax/3
sigma = 0.01
q1init = 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * np.power((pipe_model.xmid - mu) / sigma, 2)) + rho0
q1init = q1init*pipe_model.A

q2init = np.zeros(pipe_model.x.shape)
q2init = q2init[0:-1]

solq1,solq2,time = pipe_model.solve(q1init,q2init,FinalTime=FinalTime,implicit=False,stepsize=1e-6)



plt.figure()
plt.plot(solq1[-1])
plt.plot(solq2[-1])
plt.show()

animateSolution(pipe_model.xmid,time[0:-1:20],solq1[0:-1:20],gif_name='pipe_flow_simulation')

'''
plt.figure()
plt.loglog(vec,error)
plt.loglog(vec,np.divide(1,vec))
plt.grid()
plt.show()
'''