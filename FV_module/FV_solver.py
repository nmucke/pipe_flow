import numpy as np
import matplotlib.pyplot as plt
import pdb

import FV_routines
import time_integrators


class FV_solver(FV_routines.FV_1D):
    def __init__(self,xmin=0,xmax=1, num_volumes=10,
                 integrator='BDF2'):

        super(FV_solver,self).__init__(xmin=xmin,xmax=xmax,num_volumes=num_volumes)

        self.integrator = integrator
        if integrator == 'ImplicitEuler':
            self.integrator_func = self.integrator_func = time_integrators.ImplicitEuler()
        elif integrator == 'LowStorageRK':
            self.integrator_func = time_integrators.LowStorageRK()
        elif integrator == 'BDF2':
            self.integrator_func = time_integrators.BDF2()


    def solve_pde(self,q_init,t_end,step_size,rhs):
        """Solve PDE from given initial condition"""
        q_sol = [q_init]
        t_vec = [0]

        t = 0

        if self.integrator == 'BDF2':
            q_new, t_new = self.integrator_func.initial_step(time=0,
                                                      q_init=q_init,
                                                      rhs=rhs,
                                                      step_size=step_size)
            q_sol.append(q_new)
            t_vec.append(t_new)

        while t < t_end:

            if self.integrator == 'LowStorageRK':
                #h = np.dot(self.I_p, q_sol[-1][0:self.num_volumes])
                #u = q_sol[-1][-(self.num_volumes + 1):] / h

                #C = np.max(u + np.sqrt(self.g * h))
                #C = 0.001

                rhoA = np.dot(self.I_p, q_sol[-1][0:self.num_volumes])
                u = q_sol[-1][-(self.num_volumes + 1):] / rhoA

                C = self.c/np.sqrt(self.A)

                CFL = 0.5

                step_size = CFL * self.dx/C

                q_new, t_new = self.integrator_func.update_state(q_sol=q_sol,
                                                                 t_vec=t_vec,
                                                                 step_size=step_size,
                                                                 rhs=rhs)
            else:
                q_new, t_new = self.integrator_func.update_state(q_sol=q_sol,
                                                                 t_vec=t_vec,
                                                                 step_size=step_size,
                                                                 rhs=rhs)



            t = t_new

            q_sol.append(q_new)
            t_vec.append(t_new)

        return q_sol, t_vec











