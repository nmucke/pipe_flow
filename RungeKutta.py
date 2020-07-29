import numpy as np
import pdb
import scipy.linalg as scilin
import matplotlib.pyplot as plt
#import DG_routines as DG
import scipy.optimize as opt
import scipy.sparse.linalg as spla
import time as timing



class Onestepmethod(object):
    def __init__(self,f,y0,t0,te,N,jacobian,tol):
        self.f =  f
        self.y0 = y0.astype(float)
        self.t0 = t0
        self.interval = [t0, te]
        self.grid = np.linspace(t0, te, N + 2)  # N interior points
        self.h = (te - t0) / (N + 1)
        self.N = N
        self.tol = tol
        self.m = len(y0)
        self.s = len(self.b)
        self.jacobian = jacobian

    def step(self):
        t, y = [self.grid[0]], [self.y0]  # initial condition
        tim1 = t[0]
        yi = y[0]
        for ti in self.grid[1:]:
            yi = yi + self.h * self.phi(tim1, yi)
            tim1 = ti
            t.append(tim1)
            y.append(yi)

        return t, np.transpose(np.asarray(y))
    def solve(self):
        self.time, self.solution = self.step()

class RungeKutta_implicit(Onestepmethod):
    def phi(self, t0, y0):
        """
        Calculates the summation of b_j*Y_j in one step of the RungeKutta method with
        y_{n+1} = y_{n} + h*sum_{j=1}^{s} b_{j}*Y
        where j=1,2,...,s, and s is the number of stages, b the
        nodes, and Y the
        stage values of the method.
        Parameters:
        -------------
        t0 = float, current timestep
        y0 = 1 x m vector, the last solution y_n. Where m is the
        length
        of the initial condition y_0 of the IVP.
        """
        M = 10 # max number of newton iterations
        stageDer = np.array(self.s*[self.f(t0,y0)]) # initial value: Y'_0
        J = self.jacobian(t0, y0)
        stageVal = self.phi_solve(t0, y0, stageDer, J, M)

        lol = np.array([np.dot(self.b, stageVal.reshape(self.s, self.m)[:,j]) for j in range(self.m)])

        return np.array([np.dot(self.b, stageVal.reshape(self.s, self.m)[:,j]) for j in range(self.m)])


    def phi_solve(self, t0, y0, initVal, J, M):
        """
        This function solves the sm x sm system
        F(Y_i)=0
        by Newton’s method with an initial guess initVal.
        Parameters:
        -------------
        t0 = float, current timestep
        y0 = 1 x m vector, the last solution y_n. Where m is the
        length
        of the initial condition y_0 of the IVP.
        initVal = initial guess for the Newton iteration
        J = m x m matrix, the Jacobian matrix of f() evaluated in y_i
        M = maximal number of Newton iterations
        Returns:
        -------------
        The stage value Y_i
        """


        JJ = np.eye(self.s * self.m) - self.h * np.kron(self.A, J)
        luFactor = scilin.lu_factor(JJ)
        for i in range(M):
            initVal, norm_d = self.phi_newtonstep(t0, y0, initVal,luFactor)
            if norm_d < self.tol:
                #print('Newton converged in {} steps'.format(i))
                break
            elif i == M - 1:
                raise ValueError('The Newton iteration did not converge.')

        return initVal

    def phi_newtonstep(self, t0, y0, initVal, luFactor):
        """
        Takes one Newton step by solvning
        G’(Y_i)(Y^(n+1)_i-Y^(n)_i)=-G(Y_i)
        where
        G(Y_i) = Y_i - y_n - h*sum(a_{ij}*Y’_j) for j=1,...,s
        Parameters:
        -------------
        t0 = float, current timestep
        y0 = 1 x m vector, the last solution y_n. Where m is the
        length
        of the initial condition y_0 of the IVP.
        initVal = initial guess for the Newton iteration
        luFactor = (lu, piv) see documentation for linalg.lu_factor
        Returns:
        The difference Y^(n+1)_i-Y^(n)_i
        """

        d = scilin.lu_solve(luFactor, - self.F(initVal.flatten(),t0, y0))


        return initVal.flatten() + d, np.linalg.norm(d)

    def F(self, stageDer, t0, y0):
        """
        Returns the subtraction Y’_{i}-f(t_{n}+c_{i}*h, Y_{i}),
        where Y are
        the stage values, Y’ the stage derivatives and f the
        function of
        the IVP y’=f(t,y) that should be solved by the RK-method.
        Parameters:
        -------------
        stageDer = initial guess of the stage derivatives Y’
        t0 = float, current timestep
        y0 = 1 x m vector, the last solution y_n. Where m is the
        length
        of the initial condition y_0 of the IVP.
        """

        stageDer_new = np.empty((self.s, self.m))  # the i:th stageDer is on the i: th row
        for i in range(self.s):  # iterate over all stageDer
            stageVal = y0 + np.array([self.h * np.dot(self.A[i, :],stageDer.reshape(self.s, self.m)[:, j]) for j in range(self.m)])
            stageDer_new[i, :] = self.f(t0 + self.c[i] * self.h,stageVal)  # the ith stageDer is set on the ith row

        return stageDer - stageDer_new.reshape(-1)


class SDIRK(RungeKutta_implicit):

    def phi_solve(self, t0, y0, initVal, J, M):
        """
        This function solves F(Y_i)=0 by solving s systems of size m
        x m each.
        Newton’s method is used with an initial guess initVal.
        Parameters:
        -------------
        t0 = float, current timestep
        y0 = 1 x m vector, the last solution y_n. Where m is the
        length
        of the initial condition y_0 of the IVP.
        initVal = initial guess for the Newton iteration
        J = m x m matrix, the Jacobian matrix of f() evaluated in y_i
        M = maximal number of Newton iterations
        Returns:
        -------------
        The stage derivative Y’_i
        """
        """
        for i in range(self.s):  # solving the s mxm systems
            rhs = - self.F(initVal.flatten(), t0,y0)[i * self.m:(i + 1) * self.m] \
                  +np.sum([self.h * self.A[i, j] * np.dot(J, x[j]) for j in range(i)], axis=0)

            def RHS(y0):
                return - self.F(initVal.flatten(), t0,y0)[i * self.m:(i + 1) * self.m] \
                  +np.sum([self.h * self.A[i, j] * np.dot(J, x[j]) for j in range(i)], axis=0)

            initVal = opt.newton_krylov(RHS, y0)
        """
        JJ = np.eye(self.m) - self.h * self.A[0, 0] * J
        luFactor = scilin.lu_factor(JJ)
        for i in range(M):
            initVal, norm_d = self.phi_newtonstep(t0, y0, initVal, J,luFactor)
            if norm_d < self.tol:
                #print('Newton converged in {} steps'.format(i))
                break

            elif i == M - 1:
                raise ValueError('The Newton iteration did not converge.')

        return initVal

    def JacobianProduct(self,yk,ydelta):
        epsilon = 1/(np.linalg.norm(ydelta)+np.finfo(float).eps) * (np.finfo(float).eps/2)**(1/3)
        return 0.5*(self.f(self.t0,yk+epsilon*ydelta)-self.f(self.t0,yk-epsilon*ydelta))/epsilon

    def phi_newtonstep(self, t0, y0, initVal, J, luFactor):
        """
        Takes one Newton step by solvning
        G’(Y_i)(Y^(n+1)_i-Y^(n)_i)=-G(Y_i)
        where
        G(Y_i) = Y_i - haY’_i - y_n - h*sum(a_{ij}*Y’_j) for
        j=1,...,i-1
        Parameters:
        -------------
        t0 = float, current timestep
        y0 = 1 x m vector, the last solution y_n. Where m is the
        length
        of the initial condition y_0 of the IVP.
        initVal = initial guess for the Newton iteration
        luFactor = (lu, piv) see documentation for linalg.lu_factor
        Returns:
        The difference Y^(n+1)_i-Y^(n)_i
        """


        x = []
        for i in range(self.s):  # solving the s mxm systems
            '''
            rhs = - self.F(initVal.flatten(), t0,y0)[i * self.m:(i + 1) * self.m] \
                  +np.sum([self.h * self.A[i, j] * np.dot(J, x[j]) for j in range(i)], axis=0)
            d = scilin.lu_solve(luFactor, rhs)
            print(d)
            '''

            Jx = []

            rhs = - self.F(initVal.flatten(), t0, y0)[i * self.m:(i + 1) * self.m] \
                  + np.sum([self.h * self.A[i, j] *
                    self.JacobianProduct(initVal.flatten()[i * self.m:(i + 1) * self.m],x[j]) for j in range(i)], axis=0)

            MatProd = lambda v: np.dot(np.eye(self.m),v) - self.h * self.A[0, 0] * self.JacobianProduct(initVal.flatten()[i * self.m:(i + 1) * self.m],v)
            M = spla.LinearOperator((self.m, self.m), matvec = MatProd)

            d, exitCode = spla.gmres(M, rhs)


            x.append(d)


        return initVal + x, np.linalg.norm(x)



class SDIRK_tableau2s(SDIRK): # order 3
    p = (3 - np.sqrt(3))/6
    A = np.array([[p, 0], [1 - 2*p, p]])
    b = np.array([1/2, 1/2])
    c = np.array([p, 1 - p])

class SDIRK_tableau5s(SDIRK): #order 4
    A = np.array([[1/4, 0, 0, 0, 0], [1/2, 1/4, 0, 0, 0], [17/50,
    -1/25, 1/4, 0, 0],[371/1360, -137/2720, 15/544, 1/4,
    0],[25/24, -49/48, 125/16, -85/12, 1/4]])
    b = np.array([25/24, -49/48, 125/16, -85/12, 1/4])
    c = np.array([1/4, 3/4, 11/20, 1/2, 1])

class Gauss(RungeKutta_implicit): #order 6
    A=np.array([[5/36, 2/9 - np.sqrt(15)/15, 5/36 - np.sqrt(15)/30],[ 5/36 +
    np.sqrt(15)/24, 2/9, 5/36 - np.sqrt(15)/24],[ 5/36 + np.sqrt(15)/30,
    2/9 + np.sqrt(15)/15, 5/36]])
    b=[5/18,4/9,5/18]
    c=[1/2-np.sqrt(15)/10,1/2,1/2+np.sqrt(15)/10]

N = 1000
t0, te = 0, 1.
tol_newton = 1e-9
tol_sol = 1e-5
timeGrid = np.linspace(t0,te,N+2) #N interior points

def rhs(t,y):
    g = 9.81
    L = 1.

    RHS = np.array([y[1],-g/L*np.sin(y[0])])
    return RHS

def jacobian(t,y):
    g = 9.81
    L = 1.

    Jacobian = np.array([[0, 1], [-g/L*np.cos(y[0]), 0]])
    return Jacobian

tt = timing.time()
system = SDIRK_tableau2s(lambda t,y:rhs(t,y),np.array([np.pi/2.,10]),t0,te, N,lambda t,y:jacobian(t,y),tol_newton)

system.solve()
time,solution = system.time,system.solution
ttt = timing.time()

print('SDIRK2 time :' + str(ttt-tt))


plt.figure()
plt.plot(solution[0,:],solution[1,:],'-',linewidth=2.)




tt = timing.time()
system = SDIRK_tableau5s(lambda t,y:rhs(t,y),np.array([np.pi/2.,10]),t0,te, N,lambda t,y:jacobian(t,y),tol_newton)

system.solve()
time,solution = system.time,system.solution
ttt = timing.time()

print('SDIRK5 time :' + str(ttt-tt))


plt.plot(solution[0,:],solution[1,:],'-',linewidth=2.)



tt = timing.time()
system = Gauss(lambda t,y:rhs(t,y),np.array([np.pi/2.,10]),t0,te, N,lambda t,y:jacobian(t,y),tol_newton)

system.solve()
time,solution = system.time,system.solution
ttt = timing.time()

print('Gauss time :' + str(ttt-tt))


plt.plot(solution[0,:],solution[1,:],'--',linewidth=2.)

plt.show()