import numpy as np
import pdb
import DG_solver
import DG_routines


class Advection(DG_solver.DG_solver):
    def __init__(self, xmin=0,xmax=1,K=10,N=5,
                 integrator='BDF2',num_states=1,
                 params=None,**stabilizer,
                 ):

        super(Advection, self).__init__(xmin=xmin, xmax=xmax, K=K, N=N,
                                        integrator=integrator,
                                        num_states=num_states,
                                        stabilizer=stabilizer
                                        )

        self.velocity = params['velocity']
        self.inflow = params['inflow']
        self.inflow_noise = params['inflow_noise']
        self.Cv = params['Cv']
        self.xl = params['xl']
        self.xElementL = np.int(self.xl / self.xmax * self.K)

        self.BC_noise = np.random.normal(loc=0, scale=self.inflow_noise)

    def BoundaryConditions(self,time,q):
        '''Set boundary conditions'''

        qin = self.inflow + self.BC_noise
        qout = q[self.vmapO]

        return qin, qout

    def Leakage(self, time, xElementL, tl=0):

        f_l = np.zeros((self.x.shape))

        #for i in range(len(tl)):
            #if time >= tl[i, 0] and time < tl[i, 1]:

        l = np.zeros(self.N + 1)
        rl = 2 * (self.xl - self.VX[xElementL]) / self.deltax - 1
        for i in range(0, self.N + 1):
            l[i] = DG_routines.JacobiP(np.array([rl]), 0, 0, i)

        l = np.linalg.solve(np.transpose(self.V), l)

        f_l[:, xElementL] = self.Cv * l


        return f_l

    def rhs(self,time,q):

        nx = self.nx.flatten('F')

        lm = self.velocity * np.ones(q.shape)
        LFc = np.abs(np.maximum((lm[self.vmapM]), (lm[self.vmapP])))

        qFlux = self.velocity * q
        dq = q[self.vmapM] - q[self.vmapP]

        dqFlux = qFlux[self.vmapM] - qFlux[self.vmapP]
        dqFlux = nx * dqFlux / 2. - LFc / 2. * dq

        qin, qout = self.BoundaryConditions(time,q)

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

        f_l = self.Leakage(time, self.xElementL)

        rhsq = - self.rx * np.dot(self.Dr, qFlux) \
               + np.dot(self.LIFT, self.Fscale * dqFlux) \
               - self.rx * f_l

        rhsq = rhsq.flatten('F')

        return rhsq

    def solve(self,q_init,t_end,step_size):
        """Solve PDE from given initial condition"""
        q_sol = [q_init]
        t_vec = [0]

        t = 0

        if self.integrator == 'BDF2':
            q_new, t_new = self.integrator_func.initial_step(time=0,
                                                      q_init=q_init,
                                                      rhs=self.rhs,
                                                      step_size=step_size)
            q_sol.append(q_new)
            t_vec.append(t_new)

            t = t_new

        while t < t_end:
            self.BC_noise = np.random.normal(loc=0,scale=self.inflow_noise)

            if self.integrator == 'LowStorageRK':

                C = self.velocity
                CFL = 0.8
                step_size = CFL * self.dx/C

                q_new, t_new = self.integrator_func.update_state(q_sol=q_sol,
                                                                 t_vec=t_vec,
                                                                 step_size=step_size,
                                                                 rhs=self.rhs)
            else:
                q_new, t_new = self.integrator_func.update_state(q_sol=q_sol,
                                                                 t_vec=t_vec,
                                                                 step_size=step_size,
                                                                 rhs=self.rhs)


            t = t_new

            q_sol.append(q_new)
            t_vec.append(t_new)

        return q_sol, t_vec






class LinearPipeflow(DG_solver.DG_solver):
    def __init__(self, xmin=0,xmax=1,K=10,N=5,
                 integrator='BDF2',num_states=1,
                 params=None,**stabilizer,
                 ):

        super(LinearPipeflow, self).__init__(xmin=xmin, xmax=xmax, K=K, N=N,
                                        integrator=integrator,
                                        num_states=num_states,
                                        stabilizer=stabilizer
                                        )

        self.velocity = params['velocity']
        self.inflow = params['inflow']
        self.inflow_noise = params['inflow_noise']
        self.outPressure = params['outPressure']
        self.pamb = params['pamb']
        self.p0 = params['p0']
        self.rho0 = params['rho0']
        self.A = params['A']


        self.Cv = params['Cv']
        self.xl = params['xl']
        self.xElementL = np.int(self.xl / self.xmax * self.K)

        self.BC_noise = np.random.normal(loc=0, scale=self.inflow_noise)

    def BoundaryConditions(self,time,q1,q2):
        '''Set boundary conditions'''

        q1in = q1[0]
        q1out = self.outPressure

        q2in = self.inflow + self.BC_noise
        q2out = q2[-1]

        return q1in, q1out, q2in, q2out

    def Leakage(self, time, xElementL, pressure=0, rho=0):

        f_l = np.zeros((self.x.shape))

        l = np.zeros(self.N + 1)
        rl = 2 * (self.xl - self.VX[xElementL]) / self.deltax - 1
        for i in range(0, self.N + 1):
            l[i] = DG_routines.JacobiP(np.array([rl]), 0, 0, i)

        l = np.linalg.solve(np.transpose(self.V), l)

        pressureL = np.reshape(pressure, (self.Np, self.K), 'F')
        pressureL = self.EvaluateSol(np.array([self.xl]), pressureL)[0]
        rhoL = np.reshape(rho, (self.Np, self.K), 'F')
        rhoL = self.EvaluateSol(np.array([self.xl]), rhoL)[0]

        discharge_sqrt_coef = (pressureL - self.pamb) * rhoL
        f_l[:, xElementL] = self.Cv * np.sqrt(discharge_sqrt_coef) * l

        return f_l

    def rhs(self,time,q):

        nx = self.nx.flatten('F')

        q1 = q[0:int(self.Np*self.K)]
        q2 = q[-int(self.Np*self.K):]

        # Compute eigenvalue
        lm = self.u0 + self.velocity

        #Compute LFc
        LFc = np.abs(lm)

        # q1 flux
        q1Flux = self.velocity*q1 + self.K0*q2
        dq1 = q1[self.vmapM] - q1[self.vmapP]

        dq1Flux = q1Flux[self.vmapM] - q1Flux[self.vmapP]
        dq1Flux = nx * dq1Flux / 2. - LFc / 2. * dq1

        # q2 flux
        q2Flux = 1/self.rho0*q1 + self.u0*q2
        dq2 = q2[self.vmapM] - q2[self.vmapP]

        dq2Flux = q2Flux[self.vmapM] - q2Flux[self.vmapP]
        dq2Flux = nx * dq2Flux / 2. - LFc / 2. * dq2

        ### Boundary conditions ###
        q1in, q1out, q2in, q2out = self.BoundaryConditions(time,q1,q2)

        # Inflow conditions
        q1FluxIn = self.velocity*q1in + self.K0*q2in
        q2FluxIn = 1/self.rho0*q1in + self.u0*q2in
        lmIn = np.abs(lm) / 2
        nxIn = nx[self.mapI]
        dq1Flux[self.mapI] = nxIn * (q1Flux[self.vmapI] - q1FluxIn) / 2 \
                            - lmIn * (q1[self.vmapI] - q1in)
        dq2Flux[self.mapI] = nxIn * (q2Flux[self.vmapI] - q2FluxIn) / 2 \
                            - lmIn * (q2[self.vmapI] - q2in)

        # Outflow conditions
        q1FluxOut = self.velocity * q1out + self.K0 * q2out
        q2FluxOut = 1 / self.rho0 * q1out + self.u0 * q2out
        lmOut = np.abs(lm) / 2
        nxOut = nx[self.mapO]
        dq1Flux[self.mapI] = nxOut * (q1Flux[self.vmapO] - q1FluxOut) / 2 \
                             - lmOut * (q1[self.vmapO] - q1out)
        dq2Flux[self.mapI] = nxOut * (q2Flux[self.vmapO] - q2FluxOut) / 2 \
                             - lmOut * (q2[self.vmapO] - q2out)
        ######

        # Reshape flux vectors
        q1Flux = np.reshape(q1Flux, (self.Np, self.K), 'F')
        q2Flux = np.reshape(q2Flux, (self.Np, self.K), 'F')

        dq1Flux = np.reshape(dq1Flux, ((self.Nfp * self.Nfaces, self.K)), 'F')
        dq2Flux = np.reshape(dq2Flux, ((self.Nfp * self.Nfaces, self.K)), 'F')

        # Compute leakage
        f_l = self.Leakage(time, self.xElementL)

        # Compute RHS
        rhsq1 = - self.rx * np.dot(self.Dr, q1Flux) \
               + np.dot(self.LIFT, self.Fscale * dq1Flux) \
               - self.rx * f_l
        rhsq2 = - self.rx * np.dot(self.Dr, q2Flux) \
               + np.dot(self.LIFT, self.Fscale * dq2Flux) \

        rhsq1 = rhsq1.flatten('F')
        rhsq2 = rhsq2.flatten('F')

        return np.concatenate((rhsq1,rhsq2))

    def solve(self,q_init,t_end,step_size):
        """Solve PDE from given initial condition"""
        q_sol = [q_init]
        t_vec = [0]

        self.u0 = 0
        self.K0 = self.rho0*self.velocity**2

        t = 0

        if self.integrator == 'BDF2':
            q_new, t_new = self.integrator_func.initial_step(time=0,
                                                      q_init=q_init,
                                                      rhs=self.rhs,
                                                      step_size=step_size)
            q_sol.append(q_new)
            t_vec.append(t_new)

            t = t_new

        while t < t_end:
            self.BC_noise = np.random.normal(loc=0,scale=self.inflow_noise)

            if self.integrator == 'LowStorageRK':

                C = self.u0 + self.velocity
                CFL = 0.8
                step_size = CFL * self.dx/C

                q_new, t_new = self.integrator_func.update_state(q_sol=q_sol,
                                                                 t_vec=t_vec,
                                                                 step_size=step_size,
                                                                 rhs=self.rhs)
            else:
                q_new, t_new = self.integrator_func.update_state(q_sol=q_sol,
                                                                 t_vec=t_vec,
                                                                 step_size=step_size,
                                                                 rhs=self.rhs)


            t = t_new

            q_sol.append(q_new)
            t_vec.append(t_new)

        return q_sol, t_vec








