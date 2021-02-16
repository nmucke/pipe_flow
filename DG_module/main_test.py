import numpy as np
import matplotlib.pyplot as plt
import pdb

import DG_routines
import DG_solver
import DG_models

K = 20
N = 20

params = {'velocity': 0.1}
error = []
N_vec = range(1,7,1)
for N in N_vec:
    advection_model = DG_models.Advection(xmin=0,xmax=1,K=K,N=N,
                                   integrator='BDF2',
                                   params=params,
                                   stabilizer_type='None',Nc=10,s=36
                                   )

    xVec = np.reshape(advection_model.x, (N + 1) * K, 'F')

    FinalTime = .1

    mu = 0.4
    sigma = .075
    uinit_func = lambda x:1/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*((x-mu)/sigma)**2)
    uinit = uinit_func(xVec)


    u, time = advection_model.solve(uinit, t_end=FinalTime,step_size=0.001)


    true_end = uinit_func(xVec - params['velocity'] * time[-1])

    error.append(np.linalg.norm(u[-1] - true_end) / np.linalg.norm(true_end))


fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(xVec,uinit.flatten('F'),linewidth=2,label='u init')
ax1.plot(xVec,u[-1],linewidth=2,label='u approx')
ax1.plot(xVec,true_end,'--',linewidth=2,label='u true')
ax1.grid(True)
ax1.legend()

ax2.loglog(N_vec,error,'.-',linewidth=2,markersize=15)
ax2.grid(True)

plt.show()
