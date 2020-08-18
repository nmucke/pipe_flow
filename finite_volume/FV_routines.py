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


class FV_1D():
    def __init__(self,xmin,xmax,num_volumes):
        self.xmin = xmin
        self.xmax = xmax
        self.num_volumes = num_volumes

class Pipe1D(FV_1D):
    def __init__(self, xmin=0,xmax=1,num_volumes=100,c=1400,rho0=1000,p0=1e5, diameter=0.5):
        FV_1D.__init__(self, xmin=xmin,xmax=xmax,num_volumes=num_volumes)

        self.c = c
        self.rho0 = rho0
        self.p0 = p0
        self.diameter = diameter
        self.A = (diameter/2)**2 * np.pi