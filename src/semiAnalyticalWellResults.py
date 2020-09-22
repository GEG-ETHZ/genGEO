import numpy as np

from utils.fluidStateFromPT import FluidStateFromPT

class SemiAnalyticalWellResults(object):
    """docstring for SemiAnalyticalWellResults."""
    def __init__(self, N_dx, fluid):
        self.fluid = fluid
        # # create zero arrays of the size of N_dz+1 to iterate through all n+1 well segments.
        self.z_m            = np.zeros(N_dx+1)
        self.T_C_e          = np.zeros(N_dx+1)
        self.q              = np.zeros(N_dx+1)
        self.v_ms           = np.zeros(N_dx+1)
        self.delta_P_loss   = np.zeros(N_dx+1)
        self.T_C_f          = np.zeros(N_dx+1)
        self.P_Pa           = np.zeros(N_dx+1)
        self.h_Jkg          = np.zeros(N_dx+1)
        self.rho_kgm3       = np.zeros(N_dx+1)
        self.cp_JK          = np.zeros(N_dx+1)

    def finalState(self):
        return FluidStateFromPT(self.P_Pa[-1], self.T_C_f[-1], self.fluid)

    def end_P_Pa(self):
        return self.P_Pa[-1]

    def end_T_C(self):
        return self.T_C_f[-1]

    def end_h_Jkg(self):
        return self.h_Jkg[-1]

    # # TODO: get units and change name. do we need this?
    def getHeat(self):
        return -1. * np.sum(self.q)

    # # TODO: get units and change name. do we need this?
    def getPressureAlongWell(self):
        return self.P_Pa / 1.e6
