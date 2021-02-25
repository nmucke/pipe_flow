import numpy as np
import pdb
import DG_solver


class Advection(DG_solver.DG_solver):
    def __init__(self, xmin=0,xmax=1,K=10,N=5,
                 integrator='BDF2',params=None,**stabilizer,
                 ):

        super(Advection, self).__init__(xmin=xmin, xmax=xmax, K=K, N=N,
                                        integrator=integrator,
                                        stabilizer=stabilizer
                                        )

        self.velocity = params['velocity']
        self.inflow = params['inflow']
        self.inflow_noise = params['inflow_noise']


    def BoundaryConditions(self,q):
        '''Set boundary conditions'''

        qin = self.inflow + np.random.normal(loc=0,scale=self.inflow_noise)
        qout = q[self.vmapO]

        return qin, qout

    def rhs(self,time,q):

        nx = self.nx.flatten('F')

        lm = self.velocity * np.ones(q.shape)
        LFc = np.abs(np.maximum((lm[self.vmapM]), (lm[self.vmapP])))

        qFlux = self.velocity * q
        dq = q[self.vmapM] - q[self.vmapP]

        dqFlux = qFlux[self.vmapM] - qFlux[self.vmapP]
        dqFlux = nx * dqFlux / 2. - LFc / 2. * dq

        qin, qout = self.BoundaryConditions(q)

        qFluxIn = self.velocity * qin
        lmIn = np.abs(lm[self.vmapI]) / 2
        nxIn = nx[self.mapI]
        dqFlux[self.mapI] = nxIn * (qFlux[self.vmapI] - qFluxIn) / 2 \
                            - lmIn * (q[self.vmapI] - qin)

        qFluxOut = self.velocity * qout
        lmOut = np.abs(lm[self.vmapO]) / 2
        nxOut = nx[self.mapO]
        dqFlux[self.mapO] = nxOut * (qFlux[self.vmapO] - qFluxOut) / 2 \
                            - lmOut * (q[self.vmapO] - qout)

        qFlux = np.reshape(qFlux, (self.Np, self.K), 'F')

        dqFlux = np.reshape(dqFlux, ((self.Nfp * self.Nfaces, self.K)), 'F')

        rhsq = (-self.rx * np.dot(self.Dr, qFlux) + np.dot(self.LIFT,
                                                           self.Fscale * dqFlux))
        rhsq = rhsq.flatten('F')

        return rhsq

    def solve(self,q_init,t_end,step_size):
        return self.solve_pde(q_init,t_end,step_size,self.rhs)















