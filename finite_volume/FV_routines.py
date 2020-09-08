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


        self.D_p = -np.diag(np.ones(self.num_volumes)) + np.diag(np.ones(self.num_volumes-1),k=1)
        self.D_p[-1,0] = 1

        self.I_p = 0.5*np.diag(np.ones(self.num_volumes)) + 0.5*np.diag(np.ones(self.num_volumes - 1), k=-1)
        self.I_p[0, -1] = 0.5

        self.D_u = np.diag(np.ones(self.num_volumes)) - np.diag(np.ones(self.num_volumes-1), k=-1)
        self.D_u[0, -1] = -1

        self.I_u = 0.5 * np.diag(np.ones(self.num_volumes)) + 0.5 * np.diag(np.ones(self.num_volumes-1), k=1)
        self.I_u[-1, 0] = 0.5


        self.rk4a = np.array([0.0, -567301805773.0 / 1357537059087.0, -2404267990393.0 / 2016746695238.0,
                              -3550918686646.0 / 2091501179385.0, -1275806237668.0 / 842570457699.0])
        self.rk4b = np.array(
            [1432997174477.0 / 9575080441755.0, 5161836677717.0 / 13612068292357.0, 1720146321549.0 / 2090206949498.0,
             3134564353537.0 / 4481467310338.0, 2277821191437.0 / 14882151754819.0])
        self.rk4c = np.array([0.0, 1432997174477.0 / 9575080441755.0, 2526269341429.0 / 6820363962896.0,
                              2006345519317.0 / 3224310063776.0, 2802321613138.0 / 2924317926251.0])

    def minmod(self,a,b):

        sigma = []

        for i in range(len(a)):
            if np.abs(a[i]) < np.abs(b[i]) and a[i]*b[i] >0:
                sigma.append(a)
            elif np.abs(b[i]) < np.abs(a[i]) and a[i]*b[i] >0:
                sigma.append(b)
            elif a[i]*b[i]<=0:
                sigma.append(0)

        return np.asarray(sigma)

    def LF_flux(self,u_left,u_right,F_left,F_right,a):

        F = 0.5 * (F_left + F_right - a *(u_right-u_left))

        return F


class Advection1D(FV_1D):
    def __init__(self, xmin=0,xmax=1,num_volumes=100,c=1):
        FV_1D.__init__(self, xmin=xmin,xmax=xmax,num_volumes=num_volumes)

        self.c = c

    def flux(self,u):

        F = self.c*u

        return F

    def RHS(self,t,u):


        U = np.concatenate((np.array([0]), u, np.array([0])))

        a_left = np.maximum(np.abs(2*U[0:-2]),np.abs(2*U[1:-1]))
        a_right = np.maximum(np.abs(2*U[1:-1]),np.abs(2*U[2:]))
        F = self.flux(U)
        F_left = self.LF_flux(U[0:-2],U[1:-1],F[0:-2],F[1:-1],a_left)
        F_right = self.LF_flux(U[1:-1],U[2:],F[1:-1],F[2:],a_right)
        rhs = F_right-F_left
        rhs /= self.dx
        rhs = -rhs
        return rhs

    def solve(self,u0,FinalTime=10):

        t = 0
        time = [t]
        CFL = 0.5
        self.dt = CFL * self.dx/self.c

        sol = [u0]
        uold = u0
        while t < FinalTime:
            unew = uold + self.dt*self.RHS(t,uold)

            sol.append(unew)
            uold = unew
            t += self.dt
            time.append(t)

        return sol,time

    def solve_RK(self,u0,FinalTime=10):

        t = 0
        time = [t]
        CFL = 0.5
        dt = CFL * self.dx/self.c

        resu = np.zeros(u0.shape)

        sol = [u0]
        u = u0
        while t < FinalTime:
            for INTRK in range(0, 5):
                rhsu = self.RHS(t,u)

                resu = self.rk4a[INTRK] * resu + dt * rhsu

                u = u + self.rk4b[INTRK] * resu

            sol.append(u)
            t += dt
            time.append(t)

        return sol,time



class Burgers1D(FV_1D):
    def __init__(self, xmin=0,xmax=1,num_volumes=100,c=1, boundary_conditions = np.array([[0],[0]])):
        FV_1D.__init__(self, xmin=xmin,xmax=xmax,num_volumes=num_volumes)

        self.c = c
        self.boundary_conditions = boundary_conditions

    def flux(self,u):

        F = self.c*u*u

        return F

    def RHS(self,t,u):


        #inflow = np.sin(t)
        #U = np.concatenate((u[-1:],u,u[:1]))
        leftBC = lambda t: self.x[0]/(t+1)
        rightBC = lambda t: self.x[-1]/(t+1)

        left_ghost_point = np.array([leftBC(t+self.dx/(2*u[0]))])
        right_ghost_point = np.array([rightBC(t-self.dx/(2*u[-1]))])

        U = np.concatenate((left_ghost_point, u, right_ghost_point))

        a_left = np.maximum(np.abs(2*U[0:-2]),np.abs(2*U[1:-1]))
        a_right = np.maximum(np.abs(2*U[1:-1]),np.abs(2*U[2:]))
        F = self.flux(U)
        F_left = self.LF_flux(U[0:-2],U[1:-1],F[0:-2],F[1:-1],a_left)
        F_right = self.LF_flux(U[1:-1],U[2:],F[1:-1],F[2:],a_right)
        rhs = F_right-F_left
        rhs /= self.dx
        rhs = -rhs
        return rhs

    def solve(self,u0,FinalTime=10):

        t = 0
        time = [t]

        sol = [u0]
        uold = u0
        while t < FinalTime:
            CFL = 0.5
            self.dt = CFL * self.dx/np.max(2*np.abs(uold))
            unew = uold + self.dt*self.RHS(t,uold)

            sol.append(unew)
            uold = unew
            t += self.dt
            time.append(t)

        return sol,time

    def solve_RK(self,u0,FinalTime=10):

        t = 0
        time = [t]
        CFL = 0.5
        dt = CFL * self.dx/self.c

        resu = np.zeros(u0.shape)

        sol = [u0]
        u = u0
        while t < FinalTime:
            for INTRK in range(0, 5):
                rhsu = self.RHS(t,u)

                resu = self.rk4a[INTRK] * resu + dt * rhsu

                u = u + self.rk4b[INTRK] * resu

            sol.append(u)
            t += dt
            time.append(t)

        return sol,time


class LSWE(FV_1D):
    def __init__(self, xmin=0,xmax=1,num_volumes=100,d0=1,g=9.81, boundary_conditions = np.array([[0],[0]])):
        FV_1D.__init__(self, xmin=xmin,xmax=xmax,num_volumes=num_volumes)

        self.d0 = d0
        self.g = g
        self.boundary_conditions = boundary_conditions

    def flux(self,q1,q2):

        F1 = self.d0*q2
        F2 = self.g*q1

        return F1,F2

    def RHS(self,t,q1,q2):


        q1_interpolated = np.dot(self.I_p,q1)
        q1_diff = np.dot(self.D_p,q1)

        q2_interpolated = np.dot(self.I_u,q2)
        q2_diff = np.dot(self.D_u,q2)

        rhsq1 = -self.d0*q2_diff
        rhsq2 = -self.g*q1_diff

        rhsq1 = rhsq1/self.dx
        rhsq2 = rhsq2/self.dx

        return rhsq1, rhsq2

    def solve(self,u0,FinalTime=10):

        t = 0
        time = [t]

        sol = [u0]
        uold = u0
        while t < FinalTime:
            CFL = 0.5
            self.dt = CFL * self.dx/np.max(2*np.abs(uold))
            unew = uold + self.dt*self.RHS(t,uold)

            sol.append(unew)
            uold = unew
            t += self.dt
            time.append(t)

        return sol,time

    def solve_RK(self,q1init,q2init,FinalTime=10):


        resq1 = np.zeros(q1init.shape)
        resq2 = np.zeros(q2init.shape)

        solq1 = [q1init]
        solq2 = [q2init]

        q1 = q1init
        q2 = q2init

        t = 0
        time = [t]

        while t < FinalTime:
            CFL = 0.5
            dt = CFL * self.dx/np.abs(np.sqrt(self.d0)*np.sqrt(self.g))
            for INTRK in range(0, 5):

                rhsq1,rhsq2 = self.RHS(t,q1,q2)

                resq1 = self.rk4a[INTRK] * resq1 + dt * rhsq1
                resq2 = self.rk4a[INTRK] * resq2 + dt * rhsq2

                q1 = q1 + self.rk4b[INTRK] * resq1
                q2 = q2 + self.rk4b[INTRK] * resq2

            solq1.append(q1)
            solq2.append(q2)
            t = t + dt
            time.append(t)
            print(t)

        return solq1,solq2,time


class Pipe1D(FV_1D):
    def __init__(self, xmin=0,xmax=1,num_volumes=100,c=1400,rho0=1000,p0=1e5, diameter=0.5):
        FV_1D.__init__(self, xmin=xmin,xmax=xmax,num_volumes=num_volumes)

        self.c = c
        self.rho0 = rho0
        self.p0 = p0
        self.diameter = diameter
        self.A = (diameter/2)**2 * np.pi

    def RHS(self, t, q1, q2):
        q1_interpolated = np.dot(self.I_p, q1)
        q1_diff = np.dot(self.D_p, q1)

        q2_interpolated = np.dot(self.I_u, q2)
        q2_diff = np.dot(self.D_u, q2)

        

        rhsq1 = -q2_diff
        rhsq2 = -q1_diff

        rhsq1 = rhsq1 / self.dx
        rhsq2 = rhsq2 / self.dx

        return rhsq1, rhsq2

    def solve_RK(self,q1init,q2init,FinalTime=10):


        resq1 = np.zeros(q1init.shape)
        resq2 = np.zeros(q2init.shape)

        solq1 = [q1init]
        solq2 = [q2init]

        q1 = q1init
        q2 = q2init

        t = 0
        time = [t]

        while t < FinalTime:
            CFL = 0.5
            dt = CFL * self.dx/np.abs(np.sqrt(self.d0)*np.sqrt(self.g))
            for INTRK in range(0, 5):

                rhsq1,rhsq2 = self.RHS(t,q1,q2)

                resq1 = self.rk4a[INTRK] * resq1 + dt * rhsq1
                resq2 = self.rk4a[INTRK] * resq2 + dt * rhsq2

                q1 = q1 + self.rk4b[INTRK] * resq1
                q2 = q2 + self.rk4b[INTRK] * resq2

            solq1.append(q1)
            solq2.append(q2)
            t = t + dt
            time.append(t)
            print(t)

        return solq1,solq2,time