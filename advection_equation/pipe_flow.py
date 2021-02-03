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
import time_integrators as time_int

import scipy.optimize as opt
import scipy.sparse.linalg as spla
import scipy.integrate as integrate




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


    def f_leak(self,time,xElementL,tl,pressure=0,rho=0):

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

    def ComputeJacobian(self,q):

        m = 2*((self.N+1)*self.K)

        J = np.zeros((m,m))

        F = -self.PipeRHS1D(0,q)
        for col in range(m):
            pert = np.zeros(m)
            pert_jac = np.sqrt(np.finfo(float).eps) * np.maximum(np.abs(q[col]), 1)
            pert[col] = pert_jac

            qpert = q + pert
            Fpert = -self.PipeRHS1D(0,qpert)

            J[:,col] = (Fpert-F) / pert_jac

        return J

    def SteadyStateNewton(self,q):
        newton_tol = 1e-8
        MaxNewtonIter = 50

        LHS = self.ComputeJacobian(q)

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
        '''
        q1in = q1[self.vmapI]
        q1out = ((self.pOut - self.p0) / (self.c ** 2) + self.rho0) * self.A
        pin = self.c * self.c * (q1in / self.A - self.rho0) + self.p0
        pout = self.pOut
        q2in = (-q2[self.vmapI]/q1in + 2*self.uIn) * q1in
        q2out = q2[self.vmapO]
        '''
        pout = self.pOut
        q1in = q1[self.vmapI]
        q1out = ((self.pOut - self.p0) / (self.c ** 2) + self.rho0) * self.A
        pin = self.c * self.c * (q1in / self.A - self.rho0) + self.p0
        q2in = self.uIn * q1in
        q2out = q2[self.vmapO]
        '''
        q1in = q1[self.vmapO]
        q1out = q1[self.vmapI]
        q2in = q2[self.vmapO]
        q2out = q2[self.vmapI]
        pin = self.c * self.c * (q1in / self.A - self.rho0) + self.p0
        pout = self.c * self.c * (q1out / self.A - self.rho0) + self.p0
        '''
        return q1in, q1out, q2in, q2out, pin, pout


    def PipeRHS1D(self,time,q):

        q1 = q[0:int(len(q)/2)]
        q2 = q[-int(len(q)/2):]

        u = np.divide(q2, q1)

        pressure = self.c * self.c * (q1 / self.A - self.rho0) + self.p0
        self.f_l = self.f_leak(time, self.xElementL, self.tl, pressure=pressure, rho=q1 / self.A)

        cvel = self.c/np.sqrt(self.A)
        lm = np.abs(u) + cvel

        q1Flux = q2
        q2Flux = np.divide(np.power(q2, 2), q1) + pressure * self.A

        dq1 = q1[self.vmapM] - q1[self.vmapP]
        dq2 = q2[self.vmapM] - q2[self.vmapP]

        dq1Flux = q1Flux[self.vmapM] - q1Flux[self.vmapP]
        dq2Flux = q2Flux[self.vmapM] - q2Flux[self.vmapP]

        LFc = np.maximum((lm[self.vmapM]), (lm[self.vmapP]))

        dq1Flux = self.nx.flatten('F') * dq1Flux / 2. - LFc / 2. * dq1
        dq2Flux = self.nx.flatten('F') * dq2Flux / 2. - LFc / 2. * dq2

        if time >= self.tl[0,0] and time <= 200:
            self.uIn = -self.initInflow*0.9/180*time + self.initInflow*(1+0.9/180*20)
            self.pOut = -self.initOutPres*0.7/180*time + self.initOutPres*(1+0.7/180*20)

        q1in, q1out, q2in, q2out, pin, pout = self.BoundaryConditions(q1,q2,time)

        nx = self.nx.flatten('F')

        q1FluxIn = q2in
        q2FluxIn = np.divide(np.power(q2in, 2), q1in) + pin * self.A

        lmIn = lm[self.vmapI] / 2
        nxIn = nx[self.mapI]

        dq1Flux[self.mapI] = np.dot(nxIn, q1Flux[self.vmapI] - q1FluxIn) / 2 - np.dot(lmIn, q1[self.vmapI] - q1in)
        dq2Flux[self.mapI] = np.dot(nxIn, q2Flux[self.vmapI] - q2FluxIn) / 2 - np.dot(lmIn, q2[self.vmapI] - q2in)

        q1FluxOut = q2out
        q2FluxOut = np.divide(np.power(q2out, 2), q1out) + pout * self.A

        lmOut = lm[self.vmapO] / 2
        nxOut = nx[self.mapO]

        dq1Flux[self.mapO] = np.dot(nxOut, q1Flux[self.vmapO] - q1FluxOut) / 2 - np.dot(lmOut, q1[self.vmapO] - q1out)
        dq2Flux[self.mapO] = np.dot(nxOut, q2Flux[self.vmapO] - q2FluxOut) / 2 - np.dot(lmOut, q2[self.vmapO] - q2out)

        q1Flux = np.reshape(q1Flux, (self.Np, self.K), 'F')
        q2Flux = np.reshape(q2Flux, (self.Np, self.K), 'F')

        dq1Flux = np.reshape(dq1Flux, ((self.Nfp * self.Nfaces, self.K)), 'F')
        dq2Flux = np.reshape(dq2Flux, ((self.Nfp * self.Nfaces, self.K)), 'F')

        Red = q1 * self.diameter * np.abs(u) / self.A / self.mu
        f_friction = 0.25*1/((-1.8*np.log10((self.epsFriction/self.diameter/3.7)**1.11 + 6.9/Red))**2)
        friction_term = 0.5*self.diameter*np.pi*f_friction*q1/self.A*u*u
        friction_term = np.reshape(friction_term,(self.N+1,self.K),'F')

        rhsq1 = (-self.rx * np.dot(self.Dr, q1Flux) + np.dot(self.LIFT, self.Fscale * dq1Flux)) - self.rx * self.f_l
        rhsq2 = (-self.rx * np.dot(self.Dr, q2Flux) + np.dot(self.LIFT,self.Fscale * dq2Flux)) - friction_term

        rhsq1 = rhsq1.flatten('F')
        rhsq2 = rhsq2.flatten('F')

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
            '''
            for INTRK in range(0, 5):
                rhsq = self.PipeRHS1D(time, q)

                resq1 = self.rk4a[INTRK] * resq1 + dt * rhsq[0:int(m)]
                resq2 = self.rk4a[INTRK] * resq2 + dt * rhsq[-int(m):]

                q1 = q1 + self.rk4b[INTRK] * np.reshape(resq1, (self.Np, self.K), 'F')
                q2 = q2 + self.rk4b[INTRK] * np.reshape(resq2, (self.Np, self.K), 'F')

                q1 = DG_1D.SlopeLimitN(self, q1)
                q2 = DG_1D.SlopeLimitN(self, q2)

                q = np.concatenate((q1.flatten('F'), q2.flatten('F')), axis=0)
            '''
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

    def ImplicitIntegration(self, q1, q2):

        initCondition = np.concatenate((q1.flatten('F'), q2.flatten('F')), axis=0)

        system = time_int.BDF2(self.PipeRHS1D, initCondition, t0=0, te=self.FinalTime, stepsize=self.stepsize,xmin=self.xmin, xmax=self.xmax, K=self.K, N=self.N)
        t_vec, solution = system.solve()

        return solution[:,0:int(solution.shape[1] / 2)], solution[:,-int(solution.shape[1] / 2):], t_vec

    def solve(self, q1,q2, FinalTime,implicit=False,stepsize=1e-5,xl=0,tl=0,leak_type='mass',Cv=1,pamb=1e5,mu=1.,initInflow=2,initOutPres=5e5):

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


        q = np.concatenate((q1.flatten('F'), q2.flatten('F')), axis=0)
        q1,q2 = self.SteadyStateSolve(q)
        #q1 = q[0:int(len(q) / 2)]
        #q2 = q[-int(len(q) / 2):]




        benjamin_u = np.genfromtxt('u_initial.csv', delimiter=',', dtype=np.float64)
        ben_xu = benjamin_u[:, 0]
        ben_u = benjamin_u[:, 1]

        benjamin_p = np.genfromtxt('p_initial.csv', delimiter=',')
        ben_xp = benjamin_p[:, 0]
        ben_p = benjamin_p[:, 1]

        plt.figure()
        plt.plot(self.x.flatten('F'), np.divide(q2, q1), linewidth=2, label='Nikolaj')
        plt.plot(ben_xu, ben_u, linewidth=2, label='Benjamin')
        plt.xlabel('Position')
        plt.legend()
        plt.grid()
        plt.title('Velocity')
        plt.savefig('velocity_comparison')
        plt.show()

        pressure = self.c * self.c * (q1 / self.A - self.rho0) + self.p0
        plt.figure()
        plt.plot(self.x.flatten('F'), pressure, linewidth=2, label='Nikolaj')
        plt.plot(ben_xp, ben_p, linewidth=2, label='Benjamin')
        plt.xlabel('Position')
        plt.legend()
        plt.grid()
        plt.title('Pressure')
        plt.savefig('pressure_comparison')
        plt.show()



        benjamin_u = np.genfromtxt('u_initial.csv', delimiter=',',dtype=np.float64)
        ben_xu = benjamin_u[:, 0]
        ben_u = benjamin_u[:, 1]

        benjamin_p = np.genfromtxt('p_initial.csv', delimiter=',')
        ben_xp = benjamin_p[:, 0]
        ben_p = benjamin_p[:, 1]

        plt.figure()
        plt.plot(self.x.flatten('F'),np.divide(q2,q1),linewidth=2,label='Nikolaj')
        plt.plot(ben_xu,ben_u,linewidth=2,label='Benjamin')
        plt.xlabel('Position')
        plt.legend()
        plt.grid()
        plt.title('Velocity')
        plt.savefig('velocity_comparison')
        plt.show()

        pressure = self.c * self.c * (q1 / self.A - self.rho0) + self.p0
        plt.figure()
        plt.plot(self.x.flatten('F'), pressure,linewidth=2, label='Nikolaj')
        plt.plot(ben_xp, ben_p,linewidth=2, label='Benjamin')
        plt.xlabel('Position')
        plt.legend()
        plt.grid()
        plt.title('Pressure')
        plt.savefig('pressure_comparison')
        plt.show()

        q1 = self.SlopeLimitN(np.reshape(q1,(self.N+1,self.K),'F'))
        q2 = self.SlopeLimitN(np.reshape(q2,(self.N+1,self.K),'F'))

        t0 = timing.time()
        if implicit:
            solq1, solq2, tVec = self.ImplicitIntegration(q1, q2)
        else:
            solq1, solq2, tVec = Pipe1D.ExplicitIntegration(self,q1,q2,FinalTime)
        t1 = timing.time()

        print('Simulation finished \nTime: {:.2f} seconds'.format(t1-t0))

        return solq1,solq2,tVec





