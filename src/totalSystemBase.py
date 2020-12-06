
from src.fullSystemSolver import FullSystemSolverMinLCOEBrownfield
from src.fullSystemSolver import FullSystemSolverMinLCOEGreenfield
from src.fullSystemSolver import FullSystemSolverMaxPower

class TotalSystemBase(object):
    """TotalSystemBase."""

    def minimizeLCOEBrownfield(self):
        self.full_system_solver = FullSystemSolverMinLCOEBrownfield(self.full_system)
        output = self.full_system_solver.solve()
        if output.optMdot:
            return output
        else:
            raise Exception('No minimum LCOE Brownfield found!')


    def minimizeLCOEGreenfield(self):
        self.full_system_solver = FullSystemSolverMinLCOEGreenfield(self.full_system)
        output = self.full_system_solver.solve()
        if output.optMdot:
            return output
        else:
            raise Exception('No minimum LCOE Greenfield found!')


    def maximizePower(self):
        self.full_system_solver = FullSystemSolverMaxPower(self.full_system)
        output = self.full_system_solver.solve()
        if output.optMdot:
            return output
        else:
            raise Exception('No maximum Power found!')


    def solve(self, m_dot):
        return self.full_system.solve(m_dot = m_dot)
