import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy.special import gamma
import scipy.special as sci
import scipy.sparse as sps
import scipy.integrate as integrate
import time as timing
import scipy.linalg as scilin
import discontinuous_galerkin as DG

import scipy.optimize as opt
import scipy.sparse.linalg as spla
import scipy.integrate as integrate

class BDF2(DG.DG_1D):
    def __init__(self, f, u0, t0, te, stepsize=1e-5,xmin=0, xmax=1, K=10, N=2):

        DG.DG_1D.__init__(self, xmin=xmin, xmax=xmax, K=K, N=N)

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

        self.J = self.ComputeJacobian(self.sol[-1])
        J = self.J
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

        J = self.J#self.ComputeJacobian(self.sol[-1])

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
        self.Un[0:int((self.m/2))] = self.SlopeLimitN(np.reshape(self.Un[0:int((self.m/2))],(self.N+1,self.K),'F')).flatten('F')
        self.Un[-int((self.m/2)):] = self.SlopeLimitN(np.reshape(self.Un[-int((self.m/2)):],(self.N+1,self.K),'F')).flatten('F')
        self.sol.append(self.Un)
        self.time += self.deltat
        tVec.append(self.time)

        for i in range(self.Ntime - 1):
            self.Un = self.UpdateState()
            self.Un[0:int((self.m / 2))] = self.SlopeLimitN(np.reshape(self.Un[0:int((self.m / 2))], (self.N + 1, self.K), 'F')).flatten('F')
            self.Un[-int((self.m / 2)):] = self.SlopeLimitN(np.reshape(self.Un[-int((self.m / 2)):], (self.N + 1, self.K), 'F')).flatten('F')
            self.time += self.deltat
            tVec.append(self.time)
            self.sol.append(self.Un)

            if i % 2 == 0:
                print(str(int(i/(self.Ntime - 1)*100)) + '% Done' )
        return tVec, np.asarray(self.sol)
