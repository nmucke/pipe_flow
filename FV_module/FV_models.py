import numpy as np
import pdb
import FV_solver
import matplotlib.pyplot as plt


class NSWE(FV_solver.FV_solver):
    def __init__(self, xmin=0,xmax=1,num_volumes=10,
                 integrator='BDF2',params=None
                 ):

        super(NSWE, self).__init__(xmin=xmin, xmax=xmax,
                                        num_volumes=num_volumes,
                                        integrator=integrator
                                        )

        self.g = params['g']

        self.D_p[0, 0], self.D_p[0, 1] = 0, 0
        self.D_p[-1, -1] = 0

        self.I_p[0, 0], self.I_p[1, -1] = 1, 1

        self.D_u[-1, -1], self.D_u[-1, -2] = 0, 0
        self.D_u[0, 0] = 0

        self.I_u[-1, -1],self.I_u[-1, -2] = 1, 0


    def BoundaryConditions(self,q):
        '''Set boundary conditions'''

        return

    def rhs(self,time,q):


        q1 = q[0:self.num_mid]
        q2 = q[-self.num_faces:]

        q1_flux = np.dot(self.D_p,q2)
        q1_flux = -1/self.dx * q1_flux

        q2_flux = np.dot(self.D_u, 1/q1*np.dot(self.I_u,q2)**2 \
                         + 0.5 * self.g * q1**2)
        q2_flux = -1/self.dx*q2_flux
        
        rhsq = np.concatenate((q1_flux,q2_flux),axis=0)


        return rhsq

    def solve(self,q1_init,q2_init,t_end,step_size):
        q_init = np.concatenate((q1_init,q2_init))
        return self.solve_pde(q_init,t_end,step_size,self.rhs)


class Pipe(FV_solver.FV_solver):
    def __init__(self, xmin=0, xmax=1, num_volumes=10,
                 integrator='BDF2', params=None
                 ):
        super(Pipe, self).__init__(xmin=xmin, xmax=xmax,
                                   num_volumes=num_volumes,
                                   integrator=integrator
                                   )

        self.diameter = params['diameter']
        self.rho0 = params['rho0']
        self.pamb = params['pamb']
        self.c = params['c']
        self.p0 = params['p0']
        self.A = self.diameter**2/4*np.pi

        #self.D_p[0, 0], self.D_p[-1, -1] = -1, 1
        #self.D_p[-1, -1] = 0
        #self.D_p[-1, -2] = 0

        #self.I_p[0, 0], self.I_p[-1, -1] = 1, 1

        #self.D_u[-1, -1] = 0
        #self.D_u[0, 0] = 0
        #self.D_u[-1, -1], self.D_u[-1, -2] = -1, 0
        #self.D_u[0, 0], self.D_u[0, 1] = 1,1

        #self.I_u[-1, -1], self.I_u[-1, -2] = 1, 0

        #self.D_p[0, 0], self.D_p[0, 1] = 0, 0
        #self.D_p[-1, -1] = 0

        self.I_p[0, 0], self.I_p[-1, -1] = 1, 1

        self.D_u[-1, -1], self.D_u[-1, -2] = 0, 0
        self.D_u[0, 0],self.D_u[0, 1]  = 1, -1

        #self.I_u[-1, -1], self.I_u[-1, -2] = 1, 0
        pdb.set_trace()


    def BoundaryConditions(self, q):
        '''Set boundary conditions'''
        return

    def rhs(self, time, q):
        q1 = q[0:self.num_mid]
        q2 = q[-self.num_faces:]

        rhoA = q1
        u = q2/np.dot(self.I_p,rhoA)

        q1_flux = np.dot(self.D_p, np.dot(self.I_p,rhoA)*q2)
        q1_flux = -1 / self.dx * q1_flux

        pressure = self.c**2*(rhoA/self.A-self.rho0) + self.p0

        q2_flux = np.dot(self.D_u, 1 / q1 * np.dot(self.I_u, q2) ** 2 \
                         + pressure * self.A)
        q2_flux = -1 / self.dx * q2_flux

        rhsq = np.concatenate((q1_flux, q2_flux), axis=0)

        return rhsq

    def solve(self, q1_init, q2_init, t_end, step_size):
        q_init = np.concatenate((q1_init, q2_init))
        return self.solve_pde(q_init, t_end, step_size, self.rhs)


class Advection(FV_solver.FV_solver):
    def __init__(self, xmin=0,xmax=1,num_volumes=10,
                 integrator='BDF2',params=None
                 ):

        super(Advection, self).__init__(xmin=xmin, xmax=xmax,
                                        num_volumes=num_volumes,
                                        integrator=integrator
                                        )

        self.velocity = params['velocity']
        #self.inflow = params['inflow']
        #self.inflow_noise = params['inflow_noise']

        vec = np.zeros((num_volumes+1,1))
        vec[-1] = 1
        self.D_u = np.concatenate((self.D_u,vec),axis=1)
        #self.D_u[-1, -2] = 0
        self.D_u[0, 0], self.D_u[0, 1] = 0, 0

    def BoundaryConditions(self,q):
        '''Set boundary conditions'''

        qin = self.inflow + np.random.normal(loc=0,scale=self.inflow_noise)
        qout = q[self.vmapO]

        return qin, qout

    def rhs(self,time,q):

        rhsq = -1/self.dx * np.dot(self.D_u,q)

        return rhsq

    def solve(self,q_init,t_end,step_size):
        return self.solve_pde(q_init,t_end,step_size,self.rhs)















