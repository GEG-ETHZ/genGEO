import numpy as np

from src.fullSystemSolverPeakIterator import FullSystemSolverPeakIterator

class FullSystemSolverMinLCOEBrownfield(FullSystemSolverPeakIterator):
    """FullSystemSolver."""

    def __init__(self, system):
        super().__init__(system)

    def getTargetVar(self):
        output = self.full_system.gatherOutput()
        return output.capital_cost_model.LCOE_brownfield.LCOE

    def getDirection(self):
        return -1

class FullSystemSolverMinLCOEGreenfield(FullSystemSolverPeakIterator):
    """FullSystemSolver."""

    def __init__(self, system):
        super().__init__(system)

    def getTargetVar(self):
        output = self.full_system.gatherOutput()
        return output.capital_cost_model.LCOE_greenfield.LCOE

    def getDirection(self):
        return -1

class FullSystemSolverMaxPower(FullSystemSolverPeakIterator):
    """FullSystemSolver."""

    def __init__(self, system):
        super().__init__(system)

    def getTargetVar(self):
        output = self.full_system.gatherOutput()
        return output.energy_results.W_net

    def getDirection(self):
        return 1
