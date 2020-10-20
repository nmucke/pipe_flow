import numpy as np
import matplotlib.pyplot as plt
import pdb
import pipe_flow

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
poly = 'legendre'
integrator = True
for xl in [500]:
    N = 2
    K = 100

    xmin = 0.
    xmax = 5000.

    c = 1227#1400
    rho0 = 870#1000
    p0 = 5e5
    pamb = 1e5

    diameter = 0.508

    DG_model_pipe = pipe_flow.Pipe1D(xmin=xmin,xmax=xmax,K=K,N=N,c=c,rho0=rho0,p0=p0,
                              diameter=diameter,poly=poly)
    DG_model_pipe.StartUp()

    xVec = np.reshape(DG_model_pipe.x, (N + 1) * K, 'F')
    xl = 1548
    tl = np.array([[20.,1000.]])

    Cv = 4.05e-4
    mu = 1.04e-1
    initInflow = 2.
    initOutPres = 5e5

    FinalTime = 300#0.142857142857143

    mu1 = xmax/2
    sigma = 1
    #q1init =  100 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * np.power((DG_model_pipe.x - mu1) / sigma, 2)) + rho0
    #q1init *= DG_model_pipe.A
    pinit = p0*np.ones(DG_model_pipe.x.shape)
    q1init = ((pinit - DG_model_pipe.p0) / (DG_model_pipe.c ** 2) + DG_model_pipe.rho0) * DG_model_pipe.A#rho0*np.ones((N+1,K))*DG_model_pipe.A
    q2init = initInflow*np.ones((N+1,K))*q1init
    #q2init = 1 +  1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * np.power((DG_model_pipe.x - mu1) / sigma, 2))
    #q2init = q2init * q1init

    solq1,solq2, time = DG_model_pipe.solve(q1init,q2init, FinalTime=FinalTime,implicit=integrator,stepsize=2e0,
                                            xl=xl,tl=tl,leak_type='discharge',Cv=Cv,pamb=pamb,mu=mu,initInflow=initInflow,initOutPres=initOutPres)
    #%%
    rhoA = []
    u = []
    p = []
    pBar = []
    rho = []
    for i in range(len(solq1)):
        rhoA.append(solq1[i].flatten('F'))
        rho.append((solq1[i]/DG_model_pipe.A).flatten('F'))
        u.append(np.divide(solq2[i],solq1[i]).flatten('F'))
        p.append(np.reshape(c*c*(solq1[i]/DG_model_pipe.A-rho0)+p0, (N + 1) * K, 'F'))
        pBar.append(1e-5*np.reshape(c*c*(solq1[i]/DG_model_pipe.A-rho0)+p0, (N + 1) * K, 'F'))


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
    
    leaked_mass = []
    for t in range(len(time)):
        leak = DG_model_pipe.f_leak(time[t],DG_model_pipe.xElementL,tl,pressure=p[t],rho=rho[t])
        #leak = DG_model_pipe.rx*leak
        #leak_element = 0
        #for i in range(K):
        #    pdb.set_trace()
        #    leak_element += int.simps(leak[:, i], DG_model_pipe.x[:, i])
        leaked_mass.append(leak[0,DG_model_pipe.xElementL])

    pressure_leak = []
    for t in range(len(time)):
        mean_pressure = np.mean(np.reshape(p[t], (DG_model_pipe.Np, DG_model_pipe.K), 'F')[:, DG_model_pipe.xElementL])
        pressure_leak.append(mean_pressure)

    #print('Initial total mass: ' + str(initial_total_mass))
    #print('End total mass: ' + str(end_total_mass))
    #print('Mass difference: ' + str(initial_total_mass-end_total_mass))
    u_mat = np.asarray(u)
    p_mat = np.asarray(p)

    #xVec = DG_model_pipe.x.flatten('F')
    #plt.plot(time, u_mat[10,:],label='Leakage = ' + str(xl))
#plt.grid(True)
#plt.legend()
#plt.show()
plt.figure()
#plt.plot(xVec,u[0],linewidth=2,label='u')
#plt.plot(xVec,u[0],linewidth=2,label='rhoA init')
plt.plot(xVec,u[-1],linewidth=2,label='rhoA end')
#plt.plot(xVec,q1init.flatten('F'),linewidth=2,label='true init')
plt.grid(True)
plt.legend()
plt.title('Velocity at End')
plt.xlabel('Time (s)')
plt.ylabel('Velocity')
plt.show()

animateSolution(xVec,time[0:-1],p[0:-1])

'''
p_mat = np.asarray(p)

xVec_no_repeat = np.delete(xVec,range(N,(N+1)*K,N+1))

pp = []
for i in range(len(p)):
    p_new = p[i]
    #p_new =



[X,Y] = np.meshgrid(DG_model_pipe.x.flatten('F'),time)

p_mat = np.asarray(p)
plt.figure()
plt.contourf(X,Y,p_mat,levels=10)
plt.contour(X,Y,p_mat,levels=10,colors='black')
plt.xlabel('s [m]')
plt.ylabel('t [s]')
plt.title('Pressure')
plt.savefig('pressure_contour')
#plt.show()

u_mat = np.asarray(u)
plt.figure()
plt.contourf(X,Y,u_mat,levels=10)
plt.contour(X,Y,u_mat,levels=10,colors='black')
plt.xlabel('s [m]')
plt.ylabel('t [s]')
plt.title('Velocity')
plt.savefig('velocity_contour')
#plt.show()

'''
benjamin_data = np.genfromtxt('benjamin_data.csv',delimiter=',')
benjamin_time = benjamin_data[:,0]
benjamin_mass_leak = benjamin_data[:,1]

benjamin_data = np.genfromtxt('p_leak_N400_dt1.csv',delimiter=',')
benjamin_pressure_time = benjamin_data[:,0]
benjamin_pressure_leak = benjamin_data[:,1]

plt.figure()
plt.plot(time,pressure_leak,linewidth=2,label='Nikolaj')
plt.plot(benjamin_pressure_time,benjamin_pressure_leak,linewidth=2,label='Benjamin')
plt.grid(True)
plt.legend()
plt.title('Pressure at Leak')
plt.xlabel('Time (s)')
plt.ylabel('Pressure (bar)')
plt.savefig('pressure_at_leak')
#plt.show()

plt.figure()
plt.plot(time,leaked_mass,label='DG',linewidth=2)
plt.plot(benjamin_time,benjamin_mass_leak,label='Rosa',linewidth=2)
plt.grid(True)
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Mass Flow Through Leak (kg/s)')
plt.savefig('leaked_mass')
#plt.show()

plt.figure()
plt.plot(xVec,u[-1],linewidth=2,label='Nikolaj')
plt.grid(True)
plt.legend()
plt.title('Velocity at End')
plt.xlabel('Time (s)')
plt.ylabel('Velocity')
plt.show()

'''
plt.figure()
plt.plot(time,total_mass)
plt.grid(True)
plt.legend(['Total Mass'])
#plt.show()
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
'''

#%%
#animateSolution(xVec,time[0:-1],pBar[0:-1])

