import numpy as np
from utils.fluidStates import FluidState

class PorousReservoirResults(object):
    """PorousReservoirResults."""
    def __init__(self, fluid, wellspacing):
        self.fluid = fluid
        self.wellspacing = wellspacing
        self.dP = None
        self.end_P_Pa = None
        self.end_T_C = None
        self.heat = None
        self.reservoirT = None
        self.psi = None
        self.end_h_Jkg = None

    def finalState(self):
        return FluidState.getStateFromPT(self.end_P_Pa, self.end_T_C, self.fluid)

    # # TODO: do we need this?
    def getPressure(self):
        return self.end_P_Pa / 1.e6
