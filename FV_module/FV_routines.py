import numpy as np


class FV_1D():
    def __init__(self,xmin,xmax,num_volumes):
        self.xmin = xmin
        self.xmax = xmax
        self.num_volumes = num_volumes
        self.dx = (self.xmax-self.xmin)/self.num_volumes
        self.x = np.linspace(self.xmin,self.xmax,self.num_volumes+1)
        self.xmid = 0.5*(self.x[0:self.num_volumes] + self.x[1:self.num_volumes+1])

        self.num_faces = len(self.x)
        self.num_mid = len(self.xmid)


        self.D_p = -np.eye(num_volumes,num_volumes+1) \
                   + np.eye(num_volumes,num_volumes+1,1)
        self.I_p = 0.5*(np.eye(num_volumes + 1, num_volumes)
                        + np.eye(num_volumes + 1, num_volumes, -1))

        self.D_u = np.eye(num_volumes+1,num_volumes) \
                   - np.eye(num_volumes+1,num_volumes,-1)
        self.I_u = 0.5*(np.eye(num_volumes, num_volumes + 1)
                        + np.eye(num_volumes, num_volumes + 1, 1))


        #self.D_p = -np.eye(num_volumes, num_volumes) + np.eye(num_volumes, num_volumes, 1)
        #self.D_p[-1,0] = 1

        #self.I_p = 0.5 * (np.eye(num_volumes, num_volumes) + np.eye(num_volumes, num_volumes, -1))
        #self.I_p[0,-1] = 0.5

        #self.D_u = np.eye(num_volumes, num_volumes) - np.eye(num_volumes, num_volumes, -1)
        #self.D_u[0,-1] = -1

        #self.I_u = 0.5 * (np.eye(num_volumes, num_volumes) + np.eye(num_volumes, num_volumes, 1))
        #self.I_u[-1,0] = 0.5

        self.volume_size_mat = self.dx*np.eye(self.num_volumes)
        self.volume_size_mat_inv = 1/self.dx*np.eye(self.num_volumes)

