
import numpy as np
import sys

from scipy.optimize import root, minimize, newton


class FullSystemSolver(object):
    """FullSystemSolver provides a solver to determine minimum LCOE."""

    def __init__(self, system):
        self.full_system = system

    def minimizeFunctionBrownfield(self, initial_m_dot):

        if initial_m_dot >= (self.max_m_dot - 1e-8):
            return np.nan

        print(initial_m_dot)

        self.m_dots.append(initial_m_dot)
        if len(self.m_dots) > 10 and np.isclose([np.mean(self.m_dots)], [0.5], rtol = 1e-2):
            return False

        try:
            system = self.full_system.solve(initial_m_dot, self.time_years)
            output = self.full_system.gatherOutput()
            output_val = output.capital_cost_model.LCOE_brownfield.LCOE
        except ValueError as error:
            if str(error).find('TotalAnalyticSystemWater:ExceedsMaxProductionPumpPressure') > -1:
                self.max_m_dot = initial_m_dot
                output_val = np.nan
                print('NaN')
            else:
                raise(error)

        return output_val

    def minimizeLCOEBrownfield(self, time_years):

        initial_m_dot = .5
        self.m_dots = []
        self.max_m_dot = 1e4
        self.time_years = time_years

        sol = minimize(self.minimizeFunctionBrownfield, (initial_m_dot), method='Nelder-Mead',
        tol=1e-3)#, options={'maxfev': 50})#, 'maxiter': 15})

        # sol = minimize(self.minimizeFunctionBrownfield, initial_m_dot, method='Powell', tol=1e-3)

        # sol = newton(self.minimizeFunctionBrownfield, initial_m_dot, bounds=(0., 1e4))#, tol =  1., rtol = 1e-3)

        if sol.fun:
            return sol.x[0]
        else:
            return False

    def minimizeFunctionGreenfield(self, initial_m_dot):

        system = self.full_system.solve(initial_m_dot, self.time_years)

        output = self.full_system.gatherOutput()

        return output.capital_cost_model.LCOE_greenfield.LCOE

    def minimizeLCOEGreenfield(self, time_years):

        initial_m_dot = 1.
        self.time_years = time_years

        sol = minimize(self.minimizeFunctionGreenfield, (initial_m_dot), method='Nelder-Mead', tol=1e-2)

        return sol.x[0]

    def maximizeFunctionPower(self, initial_m_dot):

        system = self.full_system.solve(initial_m_dot, self.time_years)

        output = self.full_system.gatherOutput()

        return - output.energy_results.W_net

    def maximizePower(self, time_years):

        initial_m_dot = 1.
        self.time_years = time_years

        sol = minimize(self.maximizeFunctionPower, (initial_m_dot), method='Nelder-Mead', tol=1e-2)

        return sol.x[0]

    def gatherOutput(self):
        return self.full_system.gatherOutput()
