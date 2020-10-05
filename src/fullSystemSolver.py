
import numpy as np
import sys, re

from scipy.optimize import root, minimize, newton


class FullSystemSolver(object):
    """FullSystemSolver provides a solver to determine minimum LCOE."""

    def __init__(self, system):
        self.full_system = system

    def minimizeFunctionBrownfield(self, initial_m_dot):

        if len(self.m_dots) > 10 and initial_m_dot >= (self.max_m_dot - 1e-8):
            return np.nan

        if len(self.m_dots) > 6 and np.isclose([np.mean(self.m_dots)], [self.initial_m_dot], rtol = 1e-2):
            return False
        else:
            self.m_dots.append(initial_m_dot)

        print(initial_m_dot)

        try:
            system = self.full_system.solve(initial_m_dot, self.time_years)
            output = self.full_system.gatherOutput()
            output_val = output.capital_cost_model.LCOE_brownfield.LCOE

            self.test.append(np.array([initial_m_dot[0], output_val]))
            print(self.test)
        except ValueError as error:
            # regex = re.compile('Saturation pressure \[([+-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)) Pa\] corresponding to T \[([+-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)) K\] is within 1e-4 % of given p')
            # string = """Saturation pressure [7.2315e+06 Pa] corresponding to T [303.258 K] is within 1e-4 % of given p [7.2315e+06 Pa] : PropsSI("HMASS","P",7231504.074,"T",303.257922,"CO2")"""
            # re.match(regex, string)
            if str(error).find('TotalAnalyticSystemWater:ExceedsMaxProductionPumpPressure') > -1:
                print(str(error))
                self.max_m_dot = initial_m_dot
                output_val = np.nan
                print('NaN')
            elif str(error).find('FluidSystemCO2:TurbinePowerNegative') > -1 \
            or str(error).find('Saturation pressure ') > -1: # # TODO: fix this to regular expression error triggered in semiAnalyticalWell.py line: 123
                print(str(error))
                output_val = np.nan
                # output_val = np.interp(initial_m_dot, self.m_dots_test, self.lcoe_test)
                print('NaN')
            else:
                raise(error)

        return output_val

    def minimizeLCOEBrownfield(self, time_years):

        self.initial_m_dot = 10.
        self.m_dots = []
        self.max_m_dot = 1e4
        self.time_years = time_years

        # sol = minimize(self.minimizeFunctionBrownfield, (self.initial_m_dot), method='Nelder-Mead',
        # tol=1e-3)
        # print(sol.fun)
        # return sol.x[0]

        while True:
            # print('start')
            self.test = []
            sol = minimize(self.minimizeFunctionBrownfield, (self.initial_m_dot), method='Nelder-Mead',
            tol=1e-3)#, options={'maxfev': 50})#, 'maxiter': 15})
            print('out')
            # print(np.interp(1.57 * self.initial_m_dot, self.m_dots_test, self.lcoe_test))
            print(self.test)

            # sol = minimize(self.minimizeFunctionBrownfield, initial_m_dot, method='Powell', tol=1e-3)

            # sol = newton(self.minimizeFunctionBrownfield, initial_m_dot, bounds=(0., 1e4))#, tol =  1., rtol = 1e-3)
            # print('in here')
            # print(sol.fun)
            if not sol.fun:
                self.initial_m_dot = self.initial_m_dot * 0.1
                self.m_dots = []
                # print('reduce mDot')
            elif self.initial_m_dot < 1e-3:
                return False
            else:
                return sol.x[0]



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