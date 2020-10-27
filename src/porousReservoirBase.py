

class PorousReservoirBase(object):
    """PorousReservoirBase."""

    def __init__(self):
        self.params.transmissivity = lambda : self.params.permeability * self.params.reservoir_thickness
        self.params.T_reservoir = lambda : abs(self.params.depth) * abs(self.params.dT_dz) + self.params.T_surface_rock
        self.params.P_reservoir = lambda : abs(self.params.depth) * 1000. * self.params.g
        self.params.P_reservoir_max = lambda : abs(self.params.depth) * 2500. * self.params.g
