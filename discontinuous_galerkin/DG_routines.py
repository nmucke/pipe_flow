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



def JacobiP(x,alpha,beta,N):
    xp = x

    PL = np.zeros((N+1, len(xp)))


    gamma0 = 2**(alpha + beta + 1) / (alpha + beta + 1) * gamma(alpha + 1) * gamma(beta + 1) / gamma(alpha + beta + 1)
    PL[0,:] = 1.0 / np.sqrt(gamma0)
    if N == 0:
        return PL
    gamma1 = (alpha + 1) * (beta + 1) / (alpha + beta + 3) * gamma0

    PL[1,:] = ((alpha + beta + 2) * xp / 2 + (alpha - beta) / 2) / np.sqrt(gamma1)
    if N == 1:
        return PL[-1:,:]


    aold = 2 / (2 + alpha + beta) * np.sqrt((alpha + 1) * (beta + 1) / (alpha + beta + 3));
    for i in range(1,N):
        h1 = 2 * i + alpha + beta
        anew = 2 / (h1 + 2) * np.sqrt((i + 1) * (i + 1 + alpha + beta) * (i + 1 + alpha) * (i + 1 + beta) / (h1 + 1) / (h1 + 3))
        bnew = -(alpha**2-beta**2)/h1/(h1+2)
        PL[i+1,:] = 1/anew*(-aold*PL[i-1,:] + np.multiply((xp-bnew),PL[i,:]))
        aold = anew

    return PL[-1:,:]

def JacobiGQ(alpha,beta,N):
    x,w = sci.roots_jacobi(N,alpha,beta)
    return x,w

def JacobiGL(alpha,beta,N):
    x = np.zeros((N+1,1))

    if N==1:
        x[0]=-1
        x[1]=1
        x = x[:,0]

        return x
    x_int,w = JacobiGQ(alpha+1,beta+1,N-1)
    x = np.append(-1,np.append(x_int,1))
    return x

def Vandermonde1D(x,alpha,beta,N):
    V1D = np.zeros((len(x),N+1))

    for i in range(0,N+1):
        V1D[:,i] = JacobiP(x,alpha,beta,i)
    return V1D

def GradJacobiP(r,alpha,beta,N):
    dP = np.zeros((len(r),1))
    if N == 0:
        return dP
    else:
        dP[:,0] = np.sqrt(N*(N+alpha+beta+1))*JacobiP(r,alpha+1,beta+1,N-1)
    return dP

def GradVandermonde1D(r,alpha,beta,N):
    DVr = np.zeros((len(r),N+1))


    for i in range(0,N+1):

        DVr[:,i:(i+1)] = GradJacobiP(r,alpha,beta,i)
    return DVr

def Dmatrix1D(r,alpha,beta,N,V):
    Vr = GradVandermonde1D(r,alpha,beta,N)

    Dr = np.transpose(np.linalg.solve(np.transpose(V),np.transpose(Vr)))
    return Dr

def lift1D(Np,Nfaces,Nfp,V):
    Emat = np.zeros((Np,Nfaces*Nfp))
    Emat[0,0] = 1
    Emat[Np-1,1] = 1
    LIFT = np.dot(V,np.dot(np.transpose(V),Emat))
    return LIFT

def MeshGen1D(xmin,xmax,K):

    Nv = K+1

    VX = np.arange(1.,Nv+1.)

    for i in range(0,Nv):
        VX[i] = (xmax-xmin)*i/(Nv-1) + xmin

    EtoV = np.zeros((K,2))
    for k in range(0,K):
        EtoV[k,0] = k
        EtoV[k,1] = k+1

    return Nv, VX, K, EtoV

def GeometricFactors(x,Dr):
    xr = np.dot(Dr,x)
    J = xr
    rx = np.divide(1,J)

    return rx,J

def diracDelta(x):
    f = np.zeros(x.shape)
    f[np.argwhere((x<0.2e-1) & (x>-0.2e-1))] = 1
    return f

class DG_1D:
    def __init__(self, xmin=0,xmax=1,K=10,N=5):
        self.xmin = xmin
        self.xmax = xmax
        self.K = K
        self.N = N
        self.Np = N + 1

        self.NODETOL = 1e-10
        self.Nfp = 1
        self.Nfaces = 2

        self.rk4a = np.array([ 0.0,-567301805773.0/1357537059087.0,-2404267990393.0/2016746695238.0,-3550918686646.0/2091501179385.0,-1275806237668.0/842570457699.0])
        self.rk4b = np.array([ 1432997174477.0/9575080441755.0,5161836677717.0/13612068292357.0,1720146321549.0/2090206949498.0,3134564353537.0/4481467310338.0,2277821191437.0/14882151754819.0])
        self.rk4c = np.array([ 0.0,1432997174477.0/9575080441755.0,2526269341429.0/6820363962896.0,2006345519317.0/3224310063776.0,2802321613138.0/2924317926251.0])

        self.r = JacobiGL(0,0,self.N)
        self.V = Vandermonde1D(self.r,0,0,self.N)
        self.invV = np.linalg.inv(self.V)
        self.Dr = Dmatrix1D(self.r,0,0,self.N,self.V)
        self.M = np.transpose(np.linalg.solve(np.transpose(self.invV),self.invV))
        self.invM = np.linalg.inv(self.M)

        self.LIFT = lift1D(self.Np,self.Nfaces,self.Nfp,self.V)
        self.Nv, self.VX, self.K, self.EtoV = MeshGen1D(self.xmin, self.xmax, self.K)

        self.va = np.transpose(self.EtoV[:, 0])
        self.vb = np.transpose(self.EtoV[:, 1])
        self.x = np.ones((self.Np, 1)) * self.VX[self.va.astype(int)] + 0.5 * (np.reshape(self.r, (len(self.r), 1)) + 1) \
                * (self.VX[self.vb.astype(int)] - self.VX[self.va.astype(int)])
        self.deltax = np.min(np.abs(self.x[0, :] - self.x[-1, :]))

        fmask1 = np.where(np.abs(self.r + 1.) < self.NODETOL)[0]
        fmask2 = np.where(np.abs(self.r - 1.) < self.NODETOL)[0]

        self.Fmask = np.concatenate((fmask1, fmask2), axis=0)
        self.Fx = self.x[self.Fmask, :]


    def Normals1D(self):
            nx = np.zeros((self.Nfp * self.Nfaces, self.K))
            nx[0, :] = -1.0
            nx[1, :] = 1.0
            return nx

    def BuildMaps1D(self):

            nodeids = np.reshape(np.arange(0, self.K * self.Np), (self.Np, self.K), 'F')
            vmapM = np.zeros((self.Nfp, self.Nfaces, self.K))
            vmapP = np.zeros((self.Nfp, self.Nfaces, self.K))

            for k1 in range(0, self.K):
                for f1 in range(0, self.Nfaces):
                    vmapM[:, f1, k1] = nodeids[self.Fmask[f1], k1]

            for k1 in range(0, self.K):
                for f1 in range(0, self.Nfaces):
                    k2 = self.EtoE[k1, f1].astype(int)
                    f2 = self.EtoF[k1, f1].astype(int)

                    vidM = vmapM[:, f1, k1].astype(int)
                    vidP = vmapM[:, f2, k2].astype(int)

                    x1 = self.x[np.unravel_index(vidM, self.x.shape, 'F')]
                    x2 = self.x[np.unravel_index(vidP, self.x.shape, 'F')]

                    D = (x1 - x2) ** 2
                    if D < self.NODETOL:
                        vmapP[:, f1, k1] = vidP

            vmapP = vmapP.flatten('F')
            vmapM = vmapM.flatten('F')

            mapB = np.argwhere(vmapP == vmapM)
            vmapB = vmapM[mapB]

            mapI = 0
            mapO = self.K * self.Nfaces-1
            vmapI = 0
            vmapO = self.K * self.Np-1

            return vmapM.astype(int), vmapP.astype(int), vmapB.astype(int), mapB.astype(int), mapI,mapO,vmapI,vmapO

    def Connect1D(self):
        TotalFaces = self.Nfaces * self.K
        Nv = self.K + 1

        vn = [0, 1]

        SpFToV = sps.lil_matrix((TotalFaces, Nv))

        sk = 0
        for k in range(0, self.K):
            for face in range(0, self.Nfaces):
                SpFToV[sk, self.EtoV[k, vn[face]]] = 1.
                sk = sk + 1

        SpFToF = np.dot(SpFToV, np.transpose(SpFToV)) - sps.eye(TotalFaces)
        faces = np.transpose(np.nonzero(SpFToF))
        faces[:, [0, 1]] = faces[:, [1, 0]] + 1

        element1 = np.floor((faces[:, 0] - 1) / self.Nfaces)
        face1 = np.mod((faces[:, 0] - 1), self.Nfaces)
        element2 = np.floor((faces[:, 1] - 1) / self.Nfaces)
        face2 = np.mod((faces[:, 1] - 1), self.Nfaces)

        ind = np.ravel_multi_index(np.array([face1.astype(int), element1.astype(int)]), (self.Nfaces, self.K))
        EtoE = np.reshape(np.arange(0, self.K), (self.K, 1)) * np.ones((1, self.Nfaces))
        EtoE[np.unravel_index(ind, EtoE.shape, 'F')] = element2
        EtoF = np.ones((self.K, 1)) * np.reshape(np.arange(0, self.Nfaces), (1, self.Nfaces))
        EtoF[np.unravel_index(ind, EtoE.shape, 'F')] = face2
        return EtoE, EtoF

    def GeometricFactors(self):
        xr = np.dot(self.Dr, self.x)
        J = xr
        rx = np.divide(1, J)

        return rx, J

    def StartUp(self):

        self.EtoE, self.EtoF = DG_1D.Connect1D(self)

        self.vmapM, self.vmapP, self.vmapB, self.mapB,self.mapI,self.mapO,self.vmapI,self.vmapO = DG_1D.BuildMaps1D(self)

        self.nx = DG_1D.Normals1D(self)

        self.rx, self.J = DG_1D.GeometricFactors(self)

        self.alpha = 1

        self.Fscale = 1./(self.J[self.Fmask,:])

        self.a = 2*np.pi

    def minmod(self,v):

        m = v.shape[0]
        mfunc = np.zeros((v.shape[1],))
        s = np.sum(np.sign(v),0)/m

        ids = np.argwhere(np.abs(s)==1)

        if ids.shape[0]!=0:
            mfunc[ids] = s[ids] * np.min(np.abs(v[:,ids]),0)

        return mfunc

    def minmodB(self,v,M,h):

        mfunc = v[0,:]
        ids = np.argwhere(np.abs(mfunc) > M*h*h)

        if np.shape(ids)[0]>0:
            mfunc[ids[:,0]] = self.minmod(v[:,ids[:,0]])

        return mfunc

    def SlopeLimitLin(self,ul,xl,vm1,v0,vp1):

        ulimit = ul
        h = xl[self.Np-1,:]-xl[0,:]

        x0 = np.ones((self.Np,1))*(xl[0,:]+h/2)

        hN = np.ones((self.Np,1))*h

        ux = (2/hN) * np.dot(self.Dr,ul)

        ulimit = np.ones((self.Np,1))*v0 + (xl-x0)*(self.minmodB(np.stack((ux[0,:],np.divide((vp1-v0),h),np.divide((v0-vm1),h)),axis=0),M=1e-12,h=self.deltax))

        return ulimit

    def SlopeLimitN(self,u):

        uh = np.dot(self.invV,u)
        uh[1:self.Np,:] = 0
        uavg = np.dot(self.V,uh)
        v = uavg[0:1,:]

        ulimit = u
        eps0 = 1e-8

        ue1 = u[0,:]
        ue2 = u[-1:,:]

        vk = v
        vkm1 = np.concatenate((v[0,0:1],v[0,0:self.K-1]),axis=0)
        vkp1 = np.concatenate((v[0,1:self.K],v[0,(self.K-1):(self.K)]))

        ve1 = vk - DG_1D.minmod(self,np.concatenate((vk-ue1,vk-vkm1,vkp1-vk)))
        ve2 = vk + DG_1D.minmod(self,np.concatenate((ue2-vk,vk-vkm1,vkp1-vk)))

        ids = np.argwhere((np.abs(ve1-ue1)>eps0) | (np.abs(ve2-ue2)>eps0))[:,1]
        if ids.shape[0] != 0:

            uhl = np.dot(self.invV,u[:,ids])
            uhl[2:(self.Np+1),:] = 0
            ul = np.dot(self.V,uhl)

            ulimit[:,ids] = DG_1D.SlopeLimitLin(self, ul,self.x[:,ids],vkm1[ids],vk[0,ids],vkp1[ids])

        return ulimit


class BDF2(DG_1D):
    def __init__(self, f, u0, t0, te, stepsize=1e-5,xmin=0, xmax=1, K=10, N=2):

        DG_1D.__init__(self, xmin=xmin, xmax=xmax, K=K, N=N)

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

class Pipe1D(DG_1D):
    def __init__(self, xmin=0,xmax=1,K=10,N=5,c=1400,rho0=1000,p0=1e5, diameter=0.5):
        self.xmin = xmin
        self.xmax = xmax
        self.K=K
        self.N = N
        DG_1D.__init__(self, xmin=xmin,xmax=xmax,K=K,N=N)
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

    def BoundaryConditions(self,q1,q2):
        q1in = q1[self.vmapI]
        q1out = ((self.pOut - self.p0) / (self.c ** 2) + self.rho0) * self.A
        q2in = self.uIn * q1in
        q2out = q2[self.vmapO]
        pin = self.c * self.c * (q1in / self.A - self.rho0) + self.p0
        pout = self.pOut
        return q1in, q1out, q2in, q2out, pin, pout

    def PipeRHS1D(self,time,q):

        q1 = q[0:int(len(q)/2)]
        q2 = q[-int(len(q)/2):]

        u = np.divide(q2, q1)

        pressure = self.c * self.c * (q1 / self.A - self.rho0) + self.p0

        self.f_l = self.f_leak(time, self.xElementL, self.tl, pressure=pressure, rho=q1 / self.A)

        cvel = self.c
        lm = np.abs(np.divide(q2, q1)) + cvel

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

        q1in, q1out, q2in, q2out, pin, pout = self.BoundaryConditions(q1,q2)


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

        CFL = .1

        solq1 = [q1]
        solq2 = [q2]
        tVec = [time]


        resq1 = np.zeros((m))
        resq2 = np.zeros((m))

        i = 0

        q = np.concatenate((q1.flatten('F'), q2.flatten('F')), axis=0)

        while time < FinalTime:

            u = np.divide(q2, q1)
            lam = np.max(np.abs(np.concatenate((u + self.c, u - self.c))))
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
            q1_1 = DG_1D.SlopeLimitN(self, q1_1)
            q2_1 = DG_1D.SlopeLimitN(self, q2_1)
            q_1 = np.concatenate((q1_1.flatten('F'), q2_1.flatten('F')), axis=0)

            rhs = self.PipeRHS1D(time,q_1)
            rhsq1,rhsq2 = rhs[0:int(m)],rhs[-int(m):]
            rhsq1, rhsq2 = np.reshape(rhsq1, (self.Np, self.K), 'F'),np.reshape(rhsq2, (self.Np, self.K), 'F')
            q1_2 = (3*q1 + q1_1 + dt * rhsq1)/4
            q2_2 = (3*q2 + q2_1 + dt * rhsq2)/4
            q1_2 = DG_1D.SlopeLimitN(self, q1_2)
            q2_2 = DG_1D.SlopeLimitN(self, q2_2)
            q_2 = np.concatenate((q1_2.flatten('F'), q2_2.flatten('F')), axis=0)

            rhs = self.PipeRHS1D(time,q_2)
            rhsq1,rhsq2 = rhs[0:int(m)],rhs[-int(m):]
            rhsq1, rhsq2 = np.reshape(rhsq1, (self.Np, self.K), 'F'),np.reshape(rhsq2, (self.Np, self.K), 'F')
            q1 = (q1 + 2*q1_2 + 2*dt * rhsq1) / 3
            q2 = (q2 + 2*q2_2 + 2*dt * rhsq2) / 3
            q1 = DG_1D.SlopeLimitN(self, q1)
            q2 = DG_1D.SlopeLimitN(self, q2)
            q = np.concatenate((q1.flatten('F'), q2.flatten('F')), axis=0)

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

        #system = SDIRK(self.PipeRHS1DImplicit, initCondition, t0=0, te=self.FinalTime,order=1, stepsize=self.stepsize, xmin=self.xmin,xmax=self.xmax, K=self.K, N=self.N)
        system = BDF2(self.PipeRHS1D, initCondition, t0=0, te=self.FinalTime, stepsize=self.stepsize,xmin=self.xmin, xmax=self.xmax, K=self.K, N=self.N)
        t_vec, solution = system.solve()

        #solution = integrate.solve_ivp(self.PipeRHS1DImplicit, [0, self.FinalTime], initCondition,method='RK45')
        #return np.transpose(solution.y[0:int(solution.y.shape[0] / 2),:]), np.transpose(solution.y[-int(solution.y.shape[0] / 2):,:]), solution.t

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



        '''
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
        '''

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





