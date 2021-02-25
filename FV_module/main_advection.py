import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy.io import loadmat

import FV_models

params = {'velocity': 0.01,}

num_volumes = 100

Adv = FV_models.Advection(xmin=0,xmax=1,num_volumes=num_volumes,
                      params=params,integrator='ImplicitEuler')

mu = 0.5
sigma = 0.1
q1_init_func = lambda x:1/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*((x-mu)/sigma)**2)
q1_init = q1_init_func(Adv.x)

sol, time = Adv.solve(q1_init,t_end=.1,step_size=0.01)



plt.figure()
plt.plot(Adv.x,sol[0],'-',linewidth=1.5)
plt.plot(Adv.x,sol[3],'-',linewidth=1.5)
plt.show()