import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy.special import gamma
import scipy.special as sci
import scipy.sparse as sps
import scipy.integrate as integrate
import time as timing
import scipy.linalg as scilin

import scipy.optimize as opt
import scipy.sparse.linalg as spla
import scipy.integrate as integrate


class FV_1D():
    def __init__(self,xmin,xmax,num_volumes):
        self.xmin = xmin
        self.xmax = xmax
        self.num_volumes = num_volumes
        self.dx = (self.xmax-self.xmin)/self.num_volumes
        self.x = np.linspace(self.xmin,self.xmax,self.num_volumes+1)
        self.xmid = 0.5*(self.x[0:self.num_volumes] + self.x[1:self.num_volumes+1])


        self.shift_minus = np.diag(np.ones(self.num_volumes-1),k=-1)
        self.shift_plus = np.diag(np.ones(self.num_volumes-1),k=1)


class Advection1D(FV_1D):
    def __init__(self, xmin=0,xmax=1,num_volumes=100,c=1):
        FV_1D.__init__(self, xmin=xmin,xmax=xmax,num_volumes=num_volumes)

        self.c = c

    def flux(self,u):

        F = self.c*u

        return F

    def RHS(self,t,u):

        U = np.concatenate((u[-1:],u,u[:1]))



        F = 0.5*(self.flux(U[1:])+self.flux(U[0:-1])) - 0.5*self.c*(U[1:]-U[0:-1])
        rhs = F[1:] - F[:-1]
        rhs /= self.dx
        rhs = -rhs
        return rhs

    def solve(self,u0,FinalTime=10):

        t = 0
        CFL = 0.5
        dt = CFL * self.dx/self.c

        sol = [u0]
        uold = u0
        while t < FinalTime:
            unew = uold + dt*self.RHS(t,uold)

            sol.append(unew)
            uold = unew
            t += dt

        return sol




class Pipe1D(FV_1D):
    def __init__(self, xmin=0,xmax=1,num_volumes=100,c=1400,rho0=1000,p0=1e5, diameter=0.5):
        FV_1D.__init__(self, xmin=xmin,xmax=xmax,num_volumes=num_volumes)

        self.c = c
        self.rho0 = rho0
        self.p0 = p0
        self.diameter = diameter
        self.A = (diameter/2)**2 * np.pi