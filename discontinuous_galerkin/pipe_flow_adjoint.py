import numpy as np
import matplotlib.pyplot as plt
import pdb
import time as timing
import discontinuous_galerkin_global as DG
import time_integrators as time_int
import pipe_flow_global as pipe_flow
from scipy.sparse import csr_matrix


class Pipe1D_Adjoint(pipe_flow.Pipe1D):
    def __init__(self,Pipe1D=[]):
        self.Pipe1D = Pipe1D
        self.beta = 1e-6

    def measurementFunction(self,q):
        q1 = q[0:int(len(q) / 2)]
        q2 = q[-int(len(q) / 2):]

        y = np.zeros(2*self.Pipe1D.Np*self.Pipe1D.K)

        y[-1] = q2[-1]/q1[-1]

        return y

    def RHSAdjoint(self,t,adj,q,obs):

        J = self.Pipe1D.ComputeJacobian(t,q)
        y = self.measurementFunction(q)

        rhsAdj = -np.dot(np.transpose(J),adj) + self.beta*np.abs(y-obs)

        return rhsAdj


    def SolveAdjoint(self,q1,q2,obs):

        q = np.concatenate((q1,q2))
        adj = q
        rhsAdj = self.RHSAdjoint(0,adj,q,obs)

        adjoint_sol = 0

        return adjoint_sol