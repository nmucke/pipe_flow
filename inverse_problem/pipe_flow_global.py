import numpy as np
import matplotlib.pyplot as plt
import pdb
import time as timing
import discontinuous_galerkin_global as DG
import time_integrators as time_int

from scipy.sparse import csr_matrix


class Pipe1D(DG.DG_1D):
    def __init__(self, xmin=0,xmax=1,K=10,N=5,c=1400,rho0=1000,p0=1e5, diameter=0.5,poly='legendre'):
        self.xmin = xmin
        self.xmax = xmax
        self.K=K
        self.N = N
        DG.DG_1D.__init__(self, xmin=xmin,xmax=xmax,K=K,N=N,poly=poly)
        self.c = c
        self.rho0 = rho0
        self.p0 = p0
        self.diameter = diameter
        self.A = (diameter/2)**2 * np.pi
        self.L = xmax-xmin
        self.newton_tol = 1e-6
        self.MaxNewtonIter = 50
        self.reg = 1e-6

        self.alpha_BDF = np.array([1, -4/3, 1/3])
        self.beta = 2/3


    def Leakage(self,time,xElementL,tl,pressure=0,rho=0):

        f_l = np.zeros((self.x.shape))
        if self.leak_type == 'mass':
            for i in range(len(tl)):
                if time > tl[i,0] and time < tl[i,1]:
                    f_l[:, xElementL] = 1.
        elif self.leak_type == 'discharge':
            for i in range(len(tl)):
                if time >= tl[i, 0] and time < tl[i, 1]:
                    mean_pressure = np.mean(np.reshape(pressure,(self.Np,self.K),'F')[:,xElementL])
                    mean_rho = np.mean(np.reshape(rho,(self.Np,self.K),'F')[:,xElementL])
                    discharge_sqrt_coef = mean_rho*(mean_pressure-self.pamb)
                    if discharge_sqrt_coef < 0:
                        print('Discharge sqrt coefficient is negative!')
                        break
                    f_l[:, xElementL] = self.Cv*np.sqrt(discharge_sqrt_coef)
        return f_l

    def Friction(self,q1,q2,u):
        Red = q1 * self.diameter * np.abs(u) / self.A / self.mu
        f_friction = 0.25 * 1 / ((-1.8 * np.log10((self.epsFriction / self.diameter / 3.7) ** 1.11 + 6.9 / Red)) ** 2)
        friction_term = 0.5 * self.diameter * np.pi * f_friction * q1 / self.A * u * u
        return friction_term

    def ComputeJacobian(self,time,q):

        m = 2*((self.N+1)*self.K)

        J = np.zeros((m,m))

        F = self.PipeRHS1D(time,q)
        for col in range(m):
            pert = np.zeros(m)
            pert_jac = np.sqrt(np.finfo(float).eps) * np.maximum(np.abs(q[col]), 1)
            pert[col] = pert_jac

            qpert = q + pert
            Fpert = self.PipeRHS1D(time,qpert)

            J[:,col] = (Fpert-F) / pert_jac

        return J

    def SteadyStateNewton(self,q):
        newton_tol = 1e-8
        MaxNewtonIter = 50

        LHS = -self.ComputeJacobian(0,q)

        newton_error = 1e2
        iterations = 0
        q_old = q
        while newton_error > newton_tol and iterations < MaxNewtonIter:

            RHS = self.PipeRHS1D(0,q_old)

            delta_q = np.linalg.solve(LHS, RHS)

            q_old = q_old + delta_q

            newton_error = np.max(np.abs(delta_q))
            iterations = iterations + 1

        return q_old


    def SteadyStateSolve(self,q):

        q = self.SteadyStateNewton(q)

        q1 = q[0:int(len(q)/2)]
        q2 = q[-int(len(q)/2):]

        return q1,q2

    def BoundaryConditions(self,q1,q2,time):
        pout = self.pOut
        q1in = q1[self.vmapI]
        q1out = ((self.pOut - self.p0) / (self.c ** 2) + self.rho0) * self.A
        pin = self.c * self.c * (q1in / self.A - self.rho0) + self.p0
        q2in = self.uIn * q1in
        q2out = q2[self.vmapO]
        return q1in, q1out, q2in, q2out, pin, pout


    def PipeRHS1D(self,time,q):

        q1 = q[0:int(len(q)/2)]
        q2 = q[-int(len(q)/2):]

        u = np.divide(q2, q1)

        pressure = self.c * self.c * (q1 / self.A - self.rho0) + self.p0

        self.f_l = self.Leakage(time, self.xElementL, self.tl, pressure=pressure, rho=q1 / self.A)
        friction_term = self.Friction(q1,q2,u)

        cvel = self.c/np.sqrt(self.A)
        lm = np.abs(u) + cvel

        q1Flux = q2
        q2Flux = np.divide(np.power(q2, 2), q1) + pressure * self.A


        self.LFc = np.ones((self.Np*self.K))
        self.LFc[self.edgesIdx] = np.max(np.maximum((lm[self.vmapM]), (lm[self.vmapP])))

        if time >= self.tl[0,0] and time <= 200:
            self.uIn = -self.initInflow*0.9/180*time + self.initInflow*(1+0.9/180*20)
            self.pOut = -self.initOutPres*0.7/180*time + self.initOutPres*(1+0.7/180*20)

        q1in, q1out, q2in, q2out, pin, pout = self.BoundaryConditions(q1,q2,time)

        q1FluxIn = q2in
        q2FluxIn = np.divide(np.power(q2in, 2), q1in) + pin * self.A

        q1FluxOut = q2out
        q2FluxOut = np.divide(np.power(q2out, 2), q1out) + pout * self.A


        rhsq1 = -self.S_global.dot(q1Flux) - self.G_global.dot(q1Flux) - self.LFc / 2 *self.F_global.dot(q1) - self.invM_global.dot(self.f_l.flatten('F'))
        rhsq2 = -self.S_global.dot(q2Flux) - self.G_global.dot(q2Flux) - self.LFc / 2 *self.F_global.dot(q2)  - friction_term

        rhsq1 += q1FluxIn*self.e1BC - q1FluxOut*self.eNpKBC + self.LFc*(q1in*self.e1BC+q1out*self.eNpKBC)
        rhsq2 += q2FluxIn*self.e1BC - q2FluxOut*self.eNpKBC + self.LFc *(q2in*self.e1BC+q2out*self.eNpKBC)

        return np.concatenate((rhsq1,rhsq2),axis=0)

    def ExplicitIntegration(self,q1,q2,FinalTime):

        time = 0
        m = self.Np*self.K

        mindeltax = self.x[1,0]-self.x[0,0]#self.deltax

        CFL = .8

        solq1 = [q1]
        solq2 = [q2]
        tVec = [time]


        resq1 = np.zeros((m))
        resq2 = np.zeros((m))

        i = 0
        #filterMat = self.Filter1D(self.N, 0, 50)

        q = np.concatenate((q1.flatten('F'), q2.flatten('F')), axis=0)
        while time < FinalTime:

            u = np.divide(q2, q1)
            lam = np.max(np.abs(np.concatenate((u + self.c/np.sqrt(self.A), u - self.c/np.sqrt(self.A)))))
            C = np.max(lam)
            dt = CFL * mindeltax / C

            rhs = self.PipeRHS1D(time,q)
            rhsq1,rhsq2 = rhs[0:int(m)],rhs[-int(m):]
            rhsq1, rhsq2 = np.reshape(rhsq1, (self.Np, self.K), 'F'),np.reshape(rhsq2, (self.Np, self.K), 'F')
            q1_1 = q1 + dt*rhsq1
            q2_1 = q2 + dt*rhsq2
            #q1_1 = np.dot(filterMat,q1_1)
            #q2_1 = np.dot(filterMat,q2_1)
            q1_1 = self.SlopeLimitN( q1_1)
            q2_1 = self.SlopeLimitN( q2_1)
            q_1 = np.concatenate((q1_1.flatten('F'), q2_1.flatten('F')), axis=0)

            rhs = self.PipeRHS1D(time,q_1)
            rhsq1,rhsq2 = rhs[0:int(m)],rhs[-int(m):]
            rhsq1, rhsq2 = np.reshape(rhsq1, (self.Np, self.K), 'F'),np.reshape(rhsq2, (self.Np, self.K), 'F')
            q1_2 = (3*q1 + q1_1 + dt * rhsq1)/4
            q2_2 = (3*q2 + q2_1 + dt * rhsq2)/4
            q1_2 = self.SlopeLimitN(q1_2)
            q2_2 = self.SlopeLimitN(q2_2)
            #q1_2 = np.dot(filterMat,q1_2)
            #q2_2 = np.dot(filterMat,q2_2)
            q_2 = np.concatenate((q1_2.flatten('F'), q2_2.flatten('F')), axis=0)

            rhs = self.PipeRHS1D(time,q_2)
            rhsq1,rhsq2 = rhs[0:int(m)],rhs[-int(m):]
            rhsq1, rhsq2 = np.reshape(rhsq1, (self.Np, self.K), 'F'),np.reshape(rhsq2, (self.Np, self.K), 'F')
            q1 = (q1 + 2*q1_2 + 2*dt * rhsq1) / 3
            q2 = (q2 + 2*q2_2 + 2*dt * rhsq2) / 3
            q1 = self.SlopeLimitN( q1)
            q2 = self.SlopeLimitN( q2)
            #q1 = np.dot(filterMat,q1)
            #q2 = np.dot(filterMat,q2)
            q = np.concatenate((q1.flatten('F'), q2.flatten('F')), axis=0)

            #F = self.Filter1D(self.N,1,2)
            #q1 = np.dot(F,q1)
            #q2 = np.dot(F,q2)

            solq1.append(q1.flatten('F'))
            solq2.append(q2.flatten('F'))

            time = time + dt
            tVec.append(time)

            i += 1
            if i % 100 == 0:
                print(str(int(time/self.FinalTime*100)) + '% Done' )

        return solq1, solq2, tVec

    def InitialStep(self,f):

        self.J = self.ComputeJacobian(self.time,self.sol[-1])
        J = self.J
        LHS = 1/self.stepsize*np.eye(self.m) - J

        newton_error = 1e2
        iterations = 0
        U_old = self.sol[-1]
        while newton_error > self.newton_tol and iterations < self.MaxNewtonIter:
            RHS = -(1/self.stepsize*(U_old - self.sol[-1]) - f(self.time,U_old))

            delta_U = np.linalg.solve(LHS,RHS)

            U_old = U_old + delta_U

            newton_error = np.max(np.abs(delta_U))
            iterations = iterations + 1

        return U_old

    def UpdateState(self,f):

        J = self.J#self.ComputeJacobian(self.sol[-1])

        LHS = 1 / self.stepsize * np.eye(self.m) - self.beta*J

        newton_error = 1e2
        iterations = 0
        U_old = self.sol[-1]
        while newton_error > self.newton_tol and iterations < self.MaxNewtonIter:
            RHS = -(1 / self.stepsize * (self.alpha_BDF[0]*U_old + self.alpha_BDF[1]*self.sol[-1]+ self.alpha_BDF[2]*self.sol[-2]) - self.beta*f(self.time, U_old))

            delta_U = np.linalg.solve(LHS, RHS)

            U_old = U_old + delta_U

            newton_error = np.max(np.abs(delta_U))
            iterations = iterations + 1

        return U_old

    def ImplicitIntegration(self, q1, q2,FinalTime):

        q = np.concatenate((q1.flatten('F'), q2.flatten('F')), axis=0)

        self.time = 0
        self.t0 = 0
        self.m = 2*self.Np * self.K

        solq1 = [q1.flatten('F')]
        solq2 = [q2.flatten('F')]

        i = 0

        q = np.concatenate((q1.flatten('F'), q2.flatten('F')), axis=0)
        self.sol = [q]
        tVec = [self.t0]


        self.time += self.stepsize
        tVec.append(self.time)

        self.qn = self.InitialStep(self.PipeRHS1D)
        self.qn[0:int((self.m / 2))] = self.SlopeLimitN(np.reshape(self.qn[0:int((self.m / 2))], (self.N + 1, self.K), 'F')).flatten('F')
        self.qn[-int((self.m / 2)):] = self.SlopeLimitN(np.reshape(self.qn[-int((self.m / 2)):], (self.N + 1, self.K), 'F')).flatten('F')


        self.sol.append(self.qn)
        solq1.append(self.qn[0:int((self.m / 2))])
        solq2.append(self.qn[-int((self.m / 2)):])

        while self.time < FinalTime:
            self.time = self.time + self.stepsize
            tVec.append(self.time)

            self.qn = self.UpdateState(self.PipeRHS1D)
            self.qn[0:int((self.m / 2))] = self.SlopeLimitN(
                np.reshape(self.qn[0:int((self.m / 2))], (self.N + 1, self.K), 'F')).flatten('F')
            self.qn[-int((self.m / 2)):] = self.SlopeLimitN(
                np.reshape(self.qn[-int((self.m / 2)):], (self.N + 1, self.K), 'F')).flatten('F')

            self.sol.append(self.qn)

            solq1.append(self.qn[0:int((self.m / 2))])
            solq2.append(self.qn[-int((self.m / 2)):])


            i += 1
            if i % 1 == 0:
                print(str(int(self.time / self.FinalTime * 100)) + '% Done')


        return np.asarray(solq1), np.asarray(solq2), tVec

    def GlobalBoundaryConditionMatrices(self):

        I = np.eye(self.K)
        subdiag = np.eye(self.K, k=-1)
        supdiag = np.eye(self.K, k=1)


        self.e1BC = np.zeros(self.Np * self.K)
        self.e1BC[0] = 1
        self.e1BC = np.dot(self.invM_global, self.e1BC)

        self.eNpKBC = np.zeros(self.Np * self.K)
        self.eNpKBC[-1] = 1
        self.eNpKBC = np.dot(self.invM_global, self.eNpKBC)


        self.G_global = np.kron(subdiag, -0.5 * self.F) + np.kron(I, 0.5 * (self.E - self.H)) + np.kron(supdiag, 0.5 * self.G)
        self.F_global = np.kron(subdiag, -self.F) + np.kron(I, self.H + self.E) + np.kron(supdiag, -self.G)

        self.G_global[0: self.Np, 0: self.Np] += 0.5 * self.E
        self.F_global[0: self.Np, 0: self.Np] += self.E

        self.G_global[-self.Np:, -self.Np:] -= 0.5 * self.H
        self.F_global[-self.Np:, -self.Np:] += self.H


        self.G_global = np.dot(self.invM_global, self.G_global)
        self.F_global = np.dot(self.invM_global, self.F_global)


        self.invM_global = csr_matrix(self.invM_global)
        self.S_global = csr_matrix(self.S_global)
        self.G_global = csr_matrix(self.G_global)
        self.F_global = csr_matrix(self.F_global)





    def solve(self, q1,q2, FinalTime,implicit=False,stepsize=1e-5,xl=0,tl=0,
              leak_type='mass',Cv=1,pamb=1e5,mu=1.,initInflow=2,initOutPres=5e5):

        self.FinalTime = FinalTime
        self.stepsize = stepsize
        self.tl = tl
        self.xElementL = np.int(xl/self.xmax * self.K)
        self.leak_type = leak_type
        self.Cv = Cv
        self.pamb = pamb
        self.mu = mu
        self.epsFriction = 1e-8
        self.initInflow = initInflow
        self.initOutPres = initOutPres
        self.uIn = initInflow
        self.pOut = initOutPres
        self.edgesIdx = np.unique(np.concatenate((self.vmapM,self.vmapP)))

        self.GlobalBoundaryConditionMatrices()


        q = np.concatenate((q1.flatten('F'), q2.flatten('F')), axis=0)
        q1,q2 = self.SteadyStateSolve(q)


        q1 = self.SlopeLimitN(np.reshape(q1,(self.N+1,self.K),'F'))
        q2 = self.SlopeLimitN(np.reshape(q2,(self.N+1,self.K),'F'))

        t0 = timing.time()
        if implicit:
            solq1, solq2, tVec = self.ImplicitIntegration(q1, q2, FinalTime)
        else:
            solq1, solq2, tVec = Pipe1D.ExplicitIntegration(self,q1,q2,FinalTime)
        t1 = timing.time()

        print('Simulation finished \nTime: {:.2f} seconds'.format(t1-t0))

        return solq1,solq2,tVec

    def measurementFunction(self, q):
        q1 = q[0:int(len(q) / 2)]
        q2 = q[-int(len(q) / 2):]

        y = np.zeros(2 * self.Np * self.K)

        y[-1] = q2[-1] / q1[-1]

        return y

    def RHSAdjoint(self, t, adj, q, obs):

        J = self.ComputeJacobian(t, q)
        y = self.measurementFunction(q)

        rhsAdj = -np.dot(np.transpose(J), adj) + self.reg * np.abs(y - obs)

        return rhsAdj

    def SolveAdjoint(self, q1, q2, obs):

        q = np.concatenate((q1, q2))
        adj = q
        rhsAdj = self.RHSAdjoint(0, adj, q, obs)

        adjoint_sol = []

        while t > 0:


        return adjoint_sol





