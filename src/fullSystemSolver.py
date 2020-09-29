
from scipy.optimize import root, minimize, newton


class FullSystemSolver(object):
    """FullSystemSolver provides a solver to determine minimum LCOE."""

    def __init__(self, system):
        self.full_system = system

    def minimizeFunctionBrownfield(self, initial_m_dot):

        system = self.full_system.solve(initial_m_dot, self.time_years)

        output = self.full_system.gatherOutput()

        return output.capital_cost_model.LCOE_brownfield.LCOE

    def minimizeLCOEBrownfield(self, time_years):

        initial_m_dot = 1.
        self.time_years = time_years

        sol = minimize(self.minimizeFunctionBrownfield, (initial_m_dot), method='Nelder-Mead', tol=1e-3)

        return sol.x[0]

    def minimizeFunctionGreenfield(self, initial_m_dot):

        system = self.full_system.solve(initial_m_dot, self.time_years)

        output = self.full_system.gatherOutput()

        return output.capital_cost_model.LCOE_greenfield.LCOE

    def minimizeLCOEGreenfield(self, time_years):

        initial_m_dot = 1.
        self.time_years = time_years

        sol = minimize(self.minimizeFunctionGreenfield, (initial_m_dot), method='Nelder-Mead', tol=1e-3)

        return sol.x[0]

    def maximizeFunctionPower(self, initial_m_dot):

        system = self.full_system.solve(initial_m_dot, self.time_years)

        output = self.full_system.gatherOutput()

        return - output.energy_results.W_net

    def maximizePower(self, time_years):

        initial_m_dot = 1.
        self.time_years = time_years

        sol = minimize(self.maximizeFunctionPower, (initial_m_dot), method='Nelder-Mead', tol=1e-3)

        return sol.x[0]

    def gatherOutput(self):
        return self.full_system.gatherOutput()
