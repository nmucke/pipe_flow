import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy.io import loadmat

import FV_models


params = {'g': 9.82}
num_volumes = 10

NSWE = FV_models.NSWE(xmin=0,xmax=50,num_volumes=num_volumes,
                      params=params,integrator='ImplicitEuler')

q1_init = np.zeros(NSWE.xmid.shape)
q1_init[np.where(NSWE.xmid<=20)] = 3.5
q1_init[np.where(NSWE.xmid>20)] = 1.25
q2_init = np.zeros(NSWE.x.shape)

sol, time = NSWE.solve(q1_init,q2_init,t_end=2.5,step_size=0.1)
h = sol[-1][0:num_volumes]
u = sol[-1][-(num_volumes+1):]/np.dot(NSWE.I_p,sol[-1][0:num_volumes])



dambreak_data = loadmat('dambreakdata.mat')

plt.figure()
plt.plot(NSWE.xmid,h,'--',linewidth=1.5)
plt.plot(NSWE.x,u,'--',linewidth=1.5)
plt.plot(dambreak_data['x'][0],dambreak_data['h'][0],'-',linewidth=1.25)
plt.plot(dambreak_data['x'][0],dambreak_data['u'][0],'-',linewidth=1.25)
plt.show()