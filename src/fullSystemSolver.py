from utils.solver import Solver
import numpy as np

from src.fullSystemSolverPeakIterator import FullSystemSolverPeakIterator
from src.fullSystemSolverPeakIterator import SolverOptimizationType

from models.optimizationType import OptimizationType


class FullSystemSolver(FullSystemSolverPeakIterator):
    """FullSystemSolver."""

    def __init__(self, system):
        super().__init__(system)
        
    def getTargetVar(self, solveResult):
        if self.full_system.params.opt_mode == OptimizationType.MaximizePower:
            return solveResult.energy_results.W_net
        elif self.full_system.params.opt_mode == OptimizationType.MinimizeCost:
            return solveResult.capital_cost_model.LCOE_brownfield.LCOE
        else:
            raise NotImplementedError

    def getDirection(self):
        if self.full_system.params.opt_mode == OptimizationType.MaximizePower:
            return SolverOptimizationType.Maximize
        elif self.full_system.params.opt_mode == OptimizationType.MinimizeCost:
            return SolverOptimizationType.Minimize
        else:
            raise NotImplementedError
