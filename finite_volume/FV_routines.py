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


        #self.D_p = -np.eye(num_volumes,num_volumes+1) + np.eye(num_volumes,num_volumes+1,1)
        #self.I_p = 0.5*(np.eye(num_volumes + 1, num_volumes) + np.eye(num_volumes + 1, num_volumes, -1))

        #self.D_u = np.eye(num_volumes+1,num_volumes) - np.eye(num_volumes+1,num_volumes,-1)
        #self.I_u = 0.5*(np.eye(num_volumes, num_volumes + 1) + np.eye(num_volumes, num_volumes + 1, 1))

        self.D_p = -np.eye(num_volumes, num_volumes) + np.eye(num_volumes, num_volumes, 1)
        self.D_p[-1,0] = 1

        self.I_p = 0.5 * (np.eye(num_volumes, num_volumes) + np.eye(num_volumes, num_volumes, -1))
        self.I_p[0,-1] = 0.5

        self.D_u = np.eye(num_volumes, num_volumes) - np.eye(num_volumes, num_volumes, -1)
        self.D_u[0,-1] = -1

        self.I_u = 0.5 * (np.eye(num_volumes, num_volumes) + np.eye(num_volumes, num_volumes, 1))
        self.I_u[-1,0] = 0.5


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


        U = np.concatenate((np.array([u[-1]]), u, np.array([u[0]])))

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
    
    def solve_implicit(self,u0,FinalTime=10,stepsize=1e-2):

        initCondition = u0

        system = BDF2(self.RHS, initCondition, t0=0, te=FinalTime, stepsize=stepsize)
        time, sol = system.solve()

        return sol,time
    
    def solve_implicit(self,u0,FinalTime=10,stepsize=1e-2):

        initCondition = u0

        system = BDF2(self.RHS, initCondition, t0=0, te=FinalTime, stepsize=stepsize)
        time, sol = system.solve()

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

    def RHS(self,t,q1,q2):

        rhsq1 = -1/self.dx*self.d0*np.dot(self.D_p,q2)
        rhsq2 = -1/self.dx*self.g*np.dot(self.D_u,q1)

        rhsq2[0],rhsq2[-1] = 0,0

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

    def RHSimplicit(self, t, q):
        q1 = q[0:self.num_volumes]
        q2 = q[-(self.num_volumes):]

        #rho00 = q1[0]
        #rhoNplus1 = q1[-1]
        #u0 = -q2[0]/rho00
        #uNplus1 = -q2[-1]/rhoNplus1
        #q2[0] = -q2[1]
        #q2[-1] = -q2[-2]


        pressure = self.c*self.c*(q1/self.A-self.rho0) + self.p0
        #pressure0 = self.c*self.c*(rho00/self.A-self.rho0) + self.p0
        #pressureNplus1 = self.c*self.c*(rhoNplus1/self.A-self.rho0) + self.p0

        rhsq1 = -1/self.dx*np.dot(self.D_p,q2)
        rhsq2 = -1/self.dx*(np.dot(self.D_u,np.dot(self.I_u,q2*q2)/q1) + np.dot(self.D_u,pressure))

        #rhsq2[0] = -1/self.dx*(0.5*(q2[1]+q2[0])**2/q1[0]+pressure[0]   - (rho00*u0 + pressure0))
        #rhsq2[-1] = -1/self.dx*(rhoNplus1*uNplus1 + pressureNplus1 - (0.5*(q2[-1]+q2[-2])/q1[-1] + pressure[-1]) )

        rhs = np.concatenate((rhsq1,rhsq2),axis=0)

        return rhs

    def solve_explicit(self,q1init,q2init,FinalTime=10):

        resq1 = np.zeros(q1init.shape)
        resq2 = np.zeros(q2init.shape)

        solq1 = [q1init]
        solq2 = [q2init]

        q1 = q1init
        q2 = q2init

        t = 0
        time = [t]
        CFL = 0.9
        idx = 0
        while t < FinalTime:
            CFL = 0.9
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
            if idx%100 == 0:
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

        pressure = self.c*self.c*(q1/self.A-self.rho0) + self.p0

        rhsq1 = -1/self.dx*np.dot(self.D_p,q2)
        rhsq2 = -1/self.dx*(np.dot(self.D_u,np.dot(self.I_u,q2*q2)/q1) + np.dot(self.D_u,pressure)*self.A)

        #rhsq2[0],rhsq2[-1] = 0,0

        return rhsq1, rhsq2

    def RHSimplicit(self, t, q):
        q1 = q[0:self.num_volumes]
        q2 = q[-(self.num_volumes):]

        #rho00 = q1[0]
        #rhoNplus1 = q1[-1]
        #u0 = -q2[0]/rho00
        #uNplus1 = -q2[-1]/rhoNplus1
        #q2[0] = -q2[1]
        #q2[-1] = -q2[-2]

        pressure = self.c*self.c*(q1/self.A-self.rho0) + self.p0

        pressure = self.c*self.c*(q1/self.A-self.rho0) + self.p0
        #pressure0 = self.c*self.c*(rho00/self.A-self.rho0) + self.p0
        #pressureNplus1 = self.c*self.c*(rhoNplus1/self.A-self.rho0) + self.p0

        rhsq1 = -1/self.dx*np.dot(self.D_p,q2)
        rhsq2 = -1/self.dx*(np.dot(self.D_u,np.dot(self.I_u,q2*q2)/q1) + np.dot(self.D_u,pressure)*self.A)

        #rhsq2[0] = -1/self.dx*(0.5*(q2[1]+q2[0])**2/q1[0]+pressure[0]   - (rho00*u0 + pressure0))
        #rhsq2[-1] = -1/self.dx*(rhoNplus1*uNplus1 + pressureNplus1 - (0.5*(q2[-1]+q2[-2])/q1[-1] + pressure[-1]) )

        rhs = np.concatenate((rhsq1,rhsq2),axis=0)

        return rhs

    def solve_explicit(self,q1init,q2init,FinalTime=10):

        resq1 = np.zeros(q1init.shape)
        resq2 = np.zeros(q2init.shape)

        solq1 = [q1init]
        solq2 = [q2init]

        q1 = q1init
        q2 = q2init

        t = 0
        time = [t]
        CFL = 0.9
        idx = 0
        while t < FinalTime:
            u = np.dot(self.I_u,q2)/q1
            lam = np.max(np.abs(np.concatenate((u + self.c, u - self.c))))
            C = np.max(lam)
            dt = CFL * self.dx / C

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
            if idx%100 == 0:
                print(t)

            idx += 1

            idx += 1

        return solq1,solq2,time

    def solve_implicit(self,q1init,q2init,FinalTime=10):

        initCondition = np.concatenate((q1init, q2init), axis=0)

        system = BDF2(self.RHSimplicit, initCondition, t0=0, te=FinalTime, stepsize=self.stepsize)
        t_vec, solution = system.solve()

        return solution[:, 0:int(solution.shape[1] / 2)], solution[:, -int(solution.shape[1] / 2):], t_vec

    def solve(self,q1init,q2init,FinalTime=10,implicit=False,stepsize=1e2):


        if implicit:
            self.stepsize = stepsize
            solq1, solq2, time = self.solve_implicit(q1init, q2init, FinalTime=FinalTime)
        else:
            solq1,solq2,time = self.solve_explicit(q1init, q2init, FinalTime=FinalTime)

        return solq1,solq2,time


class BDF2():
    def __init__(self, f, u0, t0, te, stepsize=1e-5):

        self.f = f
        self.u0 = u0.astype(float)
        self.t0 = t0
        self.te = te
        self.deltat = stepsize
        self.Ntime = int((te-t0)/stepsize)
        self.m = len(u0)
        self.time = 0
        self.MaxNewtonIter = 50
        self.newton_tol = 1e-6

        self.alpha = np.array([1, -4/3, 1/3])
        self.beta = 2/3

    def ComputeJacobian(self,U):

        J = np.zeros((self.m,self.m))

        F = self.f(self.time,U)
        for col in range(self.m):
            pert = np.zeros(self.m)
            pert_jac = np.sqrt(np.finfo(float).eps) * np.maximum(np.abs(U[col]), 1)
            pert[col] = pert_jac

            Upert = U + pert

            Fpert = self.f(self.time,Upert)

            J[:,col] = (Fpert - F) / pert_jac

        return J

    def InitialStep(self):

        J = self.ComputeJacobian(self.sol[-1])

        LHS = 1/self.deltat*np.eye(self.m) - J

        newton_error = 1e2
        iterations = 0
        U_old = self.sol[-1]
        while newton_error > self.newton_tol and iterations < self.MaxNewtonIter:
            RHS = -(1/self.deltat*(U_old - self.sol[-1]) - self.f(self.time,U_old))

            delta_U = np.linalg.solve(LHS,RHS)

            U_old = U_old + delta_U

            newton_error = np.max(np.abs(delta_U))
            iterations = iterations + 1

        return U_old

    def UpdateState(self):

        J = self.ComputeJacobian(self.sol[-1])

        LHS = 1 / self.deltat * np.eye(self.m) - self.beta*J

        newton_error = 1e2
        iterations = 0
        U_old = self.sol[-1]
        while newton_error > self.newton_tol and iterations < self.MaxNewtonIter:
            RHS = -(1 / self.deltat * (self.alpha[0]*U_old + self.alpha[1]*self.sol[-1]+ self.alpha[2]*self.sol[-2]) - self.beta*self.f(self.time, U_old))

            delta_U = np.linalg.solve(LHS, RHS)

            U_old = U_old + delta_U

            newton_error = np.max(np.abs(delta_U))
            iterations = iterations + 1

        return U_old

    def solve(self):

        self.sol = [self.u0]
        tVec = [self.t0]
        self.time = self.t0

        self.Un = self.InitialStep()
        self.sol.append(self.Un)
        self.time += self.deltat
        tVec.append(self.time)

        for i in range(self.Ntime - 1):
            self.Un = self.UpdateState()
            self.time += self.deltat
            tVec.append(self.time)
            self.sol.append(self.Un)

            if i % 100 == 0:
                print(self.time)

        return tVec, np.asarray(self.sol)

    def solve_implicit(self,q1init,q2init,FinalTime=10):

        initCondition = np.concatenate((q1init, q2init), axis=0)

        system = BDF2(self.RHSimplicit, initCondition, t0=0, te=FinalTime, stepsize=self.stepsize)
        t_vec, solution = system.solve()

        return solution[:, 0:int(solution.shape[1] / 2)], solution[:, -int(solution.shape[1] / 2):], t_vec

    def solve(self,q1init,q2init,FinalTime=10,implicit=False,stepsize=1e2):

        rhsq1 = -1/self.dx*np.dot(self.D_p,q2)
        rhsq2 = -1/self.dx*(np.dot(self.D_u,np.dot(self.I_u,q2*q2)/q1) + np.dot(self.D_u,pressure))

        if implicit:
            self.stepsize = stepsize
            solq1, solq2, time = self.solve_implicit(q1init, q2init, FinalTime=FinalTime)
        else:
            solq1,solq2,time = self.solve_explicit(q1init, q2init, FinalTime=FinalTime)

        return solq1,solq2,time


class BDF2():
    def __init__(self, f, u0, t0, te, stepsize=1e-5):

        self.f = f
        self.u0 = u0.astype(float)
        self.t0 = t0
        self.te = te
        self.deltat = stepsize
        self.Ntime = int((te-t0)/stepsize)
        self.m = len(u0)
        self.time = 0
        self.MaxNewtonIter = 50
        self.newton_tol = 1e-6

        self.alpha = np.array([1, -4/3, 1/3])
        self.beta = 2/3

    def ComputeJacobian(self,U):

        J = np.zeros((self.m,self.m))

        F = self.f(self.time,U)
        for col in range(self.m):
            pert = np.zeros(self.m)
            pert_jac = np.sqrt(np.finfo(float).eps) * np.maximum(np.abs(U[col]), 1)
            pert[col] = pert_jac

            Upert = U + pert

            Fpert = self.f(self.time,Upert)

            J[:,col] = (Fpert - F) / pert_jac

        return J

    def InitialStep(self):

        J = self.ComputeJacobian(self.sol[-1])

        LHS = 1/self.deltat*np.eye(self.m) - J

        newton_error = 1e2
        iterations = 0
        U_old = self.sol[-1]
        while newton_error > self.newton_tol and iterations < self.MaxNewtonIter:
            RHS = -(1/self.deltat*(U_old - self.sol[-1]) - self.f(self.time,U_old))

            delta_U = np.linalg.solve(LHS,RHS)

            U_old = U_old + delta_U

            newton_error = np.max(np.abs(delta_U))
            iterations = iterations + 1

        return U_old

    def UpdateState(self):

        J = self.ComputeJacobian(self.sol[-1])

        LHS = 1 / self.deltat * np.eye(self.m) - self.beta*J

        newton_error = 1e2
        iterations = 0
        U_old = self.sol[-1]
        while newton_error > self.newton_tol and iterations < self.MaxNewtonIter:
            RHS = -(1 / self.deltat * (self.alpha[0]*U_old + self.alpha[1]*self.sol[-1]+ self.alpha[2]*self.sol[-2]) - self.beta*self.f(self.time, U_old))

            delta_U = np.linalg.solve(LHS, RHS)

            U_old = U_old + delta_U

            newton_error = np.max(np.abs(delta_U))
            iterations = iterations + 1

        return U_old

    def solve(self):

        self.sol = [self.u0]
        tVec = [self.t0]
        self.time = self.t0

        self.Un = self.InitialStep()
        self.sol.append(self.Un)
        self.time += self.deltat
        tVec.append(self.time)

        for i in range(self.Ntime - 1):
            self.Un = self.UpdateState()
            self.time += self.deltat
            tVec.append(self.time)
            self.sol.append(self.Un)

            if i % 100 == 0:
                print(self.time)

        return tVec, np.asarray(self.sol)