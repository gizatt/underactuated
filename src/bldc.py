import math
import numpy as np

from pydrake.systems.framework import VectorSystem
from pydrake.systems.analysis import Simulator
from pydrake.systems.framework import DiagramBuilder
from pydrake.systems.primitives import SignalLogger

def rotmat(theta):
    return np.array([[math.cos(theta), -math.sin(theta)],[math.sin(theta), math.cos(theta)]])

def cross2d(x1, x2):
    return x1[0] * x2[1] - x1[1] * x2[0]

# 
class BarMagnet():
    # Magnetized along the +yb axis
    xb = 1.
    yb = 1.
    zb = 1.
    M0 = 1.

    def __init__(self, x_len, y_len, z_len, M0):
        self.xb = x_len / 2.
        self.yb = y_len / 2.
        self.zb = z_len / 2.
        self.M0 = M0

    #http://www.pdi-berlin.de/fileadmin/rsync_intra4_www1/externe_Publikationen/2005/jap97rengel_publ0.pdf
    def ComputeFieldAtPoints(self, x, y, z):
        H_scaling = self.M0 / (4 * math.pi)
        H_x_sum = 0
        H_y_sum = 0
        H_z_sum = 0
        for k in range(1, 3):
            x_term = (x + self.xb*(-1)**k)
            for l in range(1, 3):
                y_term = (y + self.yb*(-1)**l)
                for m in range(1, 3):
                    z_term = (z + self.zb*(-1)**m)

                    sqrt_term = \
                        np.sqrt(x_term**2 + y_term**2 + z_term**2)

                    H_x_sum += (-1)**(k+l+m)*np.log(
                        z_term + sqrt_term)
                    H_z_sum += (-1)**(k+l+m)*np.log(
                        x_term + sqrt_term)

                    H_y_sum += (-1)**(k+l+m)*\
                        ((y_term*x_term)/(np.abs(y_term)*np.abs(x_term)))*\
                        np.arctan2(np.abs(x_term)*z_term,
                                np.abs(y_term)*sqrt_term)
        return np.array([H_scaling * H_x_sum,
                         -H_scaling * H_y_sum,
                         H_scaling * H_z_sum])



class BLDCPlant(VectorSystem):
    def __init__(self):
        VectorSystem.__init__(self,
            3,                           # Three inputs (voltage on each wire)
            2)                           # Two outputs (theta, dtheta)
        self._DeclareContinuousState(2)  # Same state as output.
        
        # Rotor rotational inertia
        self.I = 1.0
        # List of bar magnets on the rotor
        self.rotor_bar_magnets = []
        # List of coils on the stator
        self.stator_coils = []
