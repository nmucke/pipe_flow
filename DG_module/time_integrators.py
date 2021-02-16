import numpy as np
from scipy.special import gamma
import scipy.special as sci
import scipy.sparse as sps
import pdb
import DG_routines

class BDF2():
    def __init__(self, stabilizer, max_newton_iter=50, newton_tol=1e-6):

        self.stabilizer = stabilizer

        self.time = 0
        self.max_newton_iter = max_newton_iter
        self.newton_tol = newton_tol

        self.alpha = np.array([1, -4/3, 1/3])
        self.beta = 2/3

    def compute_jacobian(self,time,U,state_len,rhs):

        epsilon = np.finfo(float).eps

        J = np.zeros((state_len,state_len))

        F = rhs(time,U)
        for col in range(state_len):
            pert = np.zeros(state_len)
            pert_jac = np.sqrt(epsilon) * np.maximum(np.abs(U[col]), 1)
            pert[col] = pert_jac

            Upert = U + pert

            Fpert = rhs(time,Upert)

            J[:,col] = (Fpert - F) / pert_jac

        return J

    def initial_step(self, time, q_init, rhs, step_size):

        self.state_len =  q_init.shape[0]

        self.J = self.compute_jacobian(time,q_init,self.state_len,rhs)
        LHS = 1/step_size*np.eye(self.state_len) - self.J

        newton_error = 1e2
        iterations = 0
        q_old = q_init
        while newton_error > self.newton_tol and \
                iterations < self.max_newton_iter:
            RHS = -(1/step_size*(q_old - q_init) - rhs(time,q_old))

            delta_q = np.linalg.solve(LHS,RHS)

            q_old = q_old + delta_q

            newton_error = np.max(np.abs(delta_q))
            iterations = iterations + 1

        return self.stabilizer(q_old), time+step_size

    def update_state(self, q_sol, t_vec, step_size, rhs):

        #J = self.J#self.ComputeJacobian(self.sol[-1])

        LHS = 1 / step_size * np.eye(self.state_len) - self.beta*self.J

        newton_error = 1e2
        iterations = 0
        q_old = q_sol[-1]
        while newton_error > self.newton_tol and iterations < self.max_newton_iter:
            RHS = -(1 / step_size * (self.alpha[0] * q_old +
                                     self.alpha[1] * q_sol[-1] +
                                     self.alpha[2] * q_sol[-2]) -
                    self.beta * rhs(self.time, q_old))

            delta_q = np.linalg.solve(LHS, RHS)

            q_old = q_old + delta_q

            newton_error = np.max(np.abs(delta_q))
            iterations = iterations + 1

        return self.stabilizer(q_old), t_vec[-1]+step_size

