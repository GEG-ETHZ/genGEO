
from utils.fluidStates import FluidState

class PorousReservoirBase(object):
    """PorousReservoirBase."""

    def __init__(self):
        self.params.T_reservoir = lambda : abs(self.params.depth) * abs(self.params.dT_dz) + self.params.T_surface_rock
        self.params.P_reservoir = lambda : abs(self.params.depth) * 1000. * self.params.g
        self.params.P_reservoir_max = lambda : abs(self.params.depth) * 2500. * self.params.g
        self.params.P_system_min = lambda : FluidState.getPFromTQ(self.params.T_reservoir(), 0, self.params.working_fluid) + 1e5
