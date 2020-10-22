
from src.fullSystemSolver import FullSystemSolverMinLCOEBrownfield
from src.fullSystemSolver import FullSystemSolverMinLCOEGreenfield
from src.fullSystemSolver import FullSystemSolverMaxPower

class TotalSystemBase(object):
    """TotalSystemBase."""

    def minimizeLCOEBrownfield(self, time_years = 1.):
        self.full_system_solver = FullSystemSolverMinLCOEBrownfield(self.full_system)
        optMdot = self.full_system_solver.solve(time_years = time_years)
        if optMdot:
            output = self.full_system.gatherOutput()
            output.optMdot = optMdot
            return output
        else:
            raise Exception('No minimum LCOE Brownfield found!')


    def minimizeLCOEGreenfield(self, time_years = 1.):
        self.full_system_solver = FullSystemSolverMinLCOEGreenfield(self.full_system)
        optMdot = self.full_system_solver.solve(time_years = time_years)
        if optMdot:
            output = self.full_system.gatherOutput()
            output.optMdot = optMdot
            return output
        else:
            raise Exception('No minimum LCOE Greenfield found!')


    def maximizePower(self, time_years = 1.):
        self.full_system_solver = FullSystemSolverMaxPower(self.full_system)
        optMdot = self.full_system_solver.solve(time_years = time_years)
        if optMdot:
            output = self.full_system.gatherOutput()
            output.optMdot = optMdot
            return output
        else:
            raise Exception('No maximum Power found!')


    def solve(self, m_dot, time_years):
        self.full_system.solve(m_dot = m_dot, time_years = time_years)
        output = self.full_system.gatherOutput()
        return output
