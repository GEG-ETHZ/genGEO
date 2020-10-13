
import numpy as np
import sys, re

from scipy.optimize import root, minimize, newton


class FullSystemSolver(object):
    """FullSystemSolver provides a solver to determine minimum LCOE."""

    def __init__(self, system):
        self.full_system = system

    def minimizeFunctionBrownfield(self, initial_m_dot):

        if len(self.m_dots) > 10 and initial_m_dot >= (self.max_m_dot - 1e-8):
            print('NaN')
            return np.nan

        if len(self.m_dots) > 6 and np.isclose([np.mean(self.m_dots)], [self.initial_m_dot], rtol = 1e-2):
            print('STOP')
            return False
        else:
            self.m_dots.append(initial_m_dot)

        print('main val', initial_m_dot)

        try:
            system = self.full_system.solve(initial_m_dot, self.time_years)
            output = self.full_system.gatherOutput()
            output_val = output.capital_cost_model.LCOE_brownfield.LCOE

            if not np.isnan(output_val):
                self.test = np.array([initial_m_dot[0], output_val[0]])

        except Exception as error:
            print(str(error))
            # regex = re.compile('Saturation pressure \[([+-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)) Pa\] corresponding to T \[([+-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)) K\] is within 1e-4 % of given p')
            # string = """Saturation pressure [7.2315e+06 Pa] corresponding to T [303.258 K] is within 1e-4 % of given p [7.2315e+06 Pa] : PropsSI("HMASS","P",7231504.074,"T",303.257922,"CO2")"""
            # re.match(regex, string)
            if str(error).find('TotalAnalyticSystemWater:ExceedsMaxProductionPumpPressure') > -1:
                # print(str(error))
                self.max_m_dot = initial_m_dot
                output_val = np.nan
                print('NaN')
            elif str(error).find('FluidSystemCO2:TurbinePowerNegative') > -1 \
            or str(error).find('Saturation pressure ') > -1: # # TODO: fix this to regular expression error triggered in semiAnalyticalWell.py line: 123
                # print(str(error))
                n = 2
                while n <= 3:
                    try:
                        diff = initial_m_dot[0] - self.test[0]
                        tmp_m_dot = self.test[0] + n * diff
                        # print(tmp_m_dot)
                        system = self.full_system.solve(tmp_m_dot, self.time_years)
                        output = self.full_system.gatherOutput()
                        output_val_tmp = output.capital_cost_model.LCOE_brownfield.LCOE
                        output_val = np.interp(initial_m_dot, [self.test[0], tmp_m_dot], [self.test[1], output_val_tmp])[0]
                        print('tmp val ', output_val_tmp)
                        # print(output_val)
                        n = 11
                    except:
                        n += 1
                        output_val = np.nan
                        print('NaN')
            else:
                raise(error)

        # print(output_val)
        return output_val

    def minimizeLCOEBrownfield(self, time_years):

        self.initial_m_dot = 20.
        self.m_dots = []
        self.max_m_dot = 1e4
        self.time_years = time_years

        # sol = minimize(self.minimizeFunctionBrownfield, (self.initial_m_dot), method='Nelder-Mead',
        # tol=1e-3)
        # print(sol.fun)
        # return sol.x[0]

        while True:
            # sol = minimize(self.minimizeFunctionBrownfield, (self.initial_m_dot), method='TNC',
            # tol=1e-3, options={'eps': 1e0, 'scale': None, 'offset': None, 'mesg_num': None,
            # 'maxCGit': - 1, 'maxiter': None, 'eta': - 1, 'stepmx': 0, 'accuracy': 0,
            # 'minfev': 0, 'ftol': - 1, 'xtol': - 1, 'gtol': - 1, 'rescale': - 1, 'disp': False,
            # 'finite_diff_rel_step': None, 'maxfun': None})

            sol = minimize(self.minimizeFunctionBrownfield, (self.initial_m_dot),
            method='Nelder-Mead', tol=1e-3)#, options={'maxfev': 50})#, 'maxiter': 15})

            # sol = minimize(self.minimizeFunctionBrownfield, initial_m_dot, method='Powell', tol=1e-3)

            # sol = newton(self.minimizeFunctionBrownfield, initial_m_dot, bounds=(0., 1e4))#, tol =  1., rtol = 1e-3)
            # print('in here')
            # print(sol.fun)
            if not sol.fun:
                self.initial_m_dot = self.initial_m_dot * 0.1
                self.m_dots = []
                # print('reduce mDot')
                if self.initial_m_dot < 1e-3:
                    # print('STOP')
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
