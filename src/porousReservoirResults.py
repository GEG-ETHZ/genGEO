import numpy as np

from utils.fluidStateFromPT import FluidStateFromPT

class PorousReservoirResults(object):
    """PorousReservoirResults."""
    def __init__(self):
        self.dP = None
        self.heat = None
        self.psi = None
        self.state = None
