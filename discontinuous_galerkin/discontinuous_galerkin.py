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
    def __init__(self, xmin=0,xmax=1,K=10,N=5,poly='legendre'):
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

        if poly == 'legendre':
            alpha = 0
            beta = 0
        elif poly == 'chebyshev':
            alpha = -0.5
            beta = -0.5

        self.r = JacobiGL(alpha, beta, self.N)
        self.Nv, self.VX, self.K, self.EtoV = MeshGen1D(self.xmin, self.xmax, self.K)

        self.va = np.transpose(self.EtoV[:, 0])
        self.vb = np.transpose(self.EtoV[:, 1])
        self.x = np.ones((self.Np, 1)) * self.VX[self.va.astype(int)] + 0.5 * (np.reshape(self.r, (len(self.r), 1)) + 1) \
                * (self.VX[self.vb.astype(int)] - self.VX[self.va.astype(int)])
        self.deltax = np.min(np.abs(self.x[0, :] - self.x[-1, :]))

        self.r = JacobiGL(alpha,beta,self.N)
        self.V = Vandermonde1D(self.r,alpha,beta,self.N)
        self.invV = np.linalg.inv(self.V)
        self.Dr = Dmatrix1D(self.r,alpha,beta,self.N,self.V)
        #self.M = np.transpose(np.linalg.solve(np.transpose(self.invV),self.invV))
        #self.invM = np.linalg.inv(self.M)
        self.invM = np.dot(self.V, np.transpose(self.V))
        self.M = np.linalg.inv(self.invM)
        self.invMk = 2 / self.deltax * self.invM

        self.LIFT = lift1D(self.Np,self.Nfaces,self.Nfp,self.V)


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

        #ulimit = np.ones((self.Np,1))*v0 + (xl-x0)*(self.minmodB(np.stack((ux[0,:],np.divide((vp1-v0),h),np.divide((v0-vm1),h)),axis=0),M=1e-12,h=self.deltax))
        ulimit = np.ones((self.Np,1))*v0 + (xl-x0)*(self.minmod(np.stack((ux[0,:],np.divide((vp1-v0),h),np.divide((v0-vm1),h)),axis=0)))

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

    def Filter1D(self,N,Nc,s):


        filterdiag = np.ones((N+1))

        alpha = -np.log(np.finfo(float).eps)


        for i in range(Nc,N):
            #filterdiag[i+1] = np.exp(-alpha*((i-Nc)/(N-Nc))**s)
            filterdiag[i+1] = np.exp(-alpha*((i-1)/N)**s)


        F = np.dot(self.V,np.dot(np.diag(filterdiag),self.invV))
        return F

    def FindElement(self,x):

        diff = x-self.VX
        element = np.argwhere(diff >= 0)[-1,0]

        return element

    def EvaluateSol(self,x,sol_nodal):

        sol_modal = np.dot(self.invV, sol_nodal)

        if x == np.any(self.VX):

            i_interface = np.argwhere(x == self.VX)

            sol_xVec = []
            for i in range(len(x)):
                if x[i] == self.xmin:
                    sol_x = sol_nodal[0,0]
                elif x[i] == self.xmax:
                    sol_x = sol_nodal[-1, -1]
                elif i == i_interface and i != 0 and i != len(x)-1:
                    sol_x = 0.5*(sol_nodal[i-1,-1]+sol_nodal[i,0])
                else:
                    element = self.FindElement(x[i])

                    x_ref = 2 * (x[i] - self.VX[element]) / self.deltax - 1

                    sol_x = 0
                    for j in range(self.Np):
                        P = JacobiP(np.array([x_ref]), 0, 0, j)
                        sol_x += sol_modal[j, element] * P[0, 0]

                sol_xVec.append(sol_x)

        else:

            sol_xVec = []
            for i in range(len(x)):
                if x[i] == self.xmin:
                    sol_x = sol_nodal[0, 0]
                elif x[i] == self.xmax:
                    sol_x = sol_nodal[-1, -1]
                else:
                    element = self.FindElement(x[i])

                    x_ref = 2*(x[i]-self.VX[element])/self.deltax-1

                    sol_x = 0
                    for j in range(self.Np):
                        P = JacobiP(np.array([x_ref]), 0, 0, j)
                        sol_x += sol_modal[j,element]*P[0,0]

                sol_xVec.append(sol_x)

        return sol_xVec
