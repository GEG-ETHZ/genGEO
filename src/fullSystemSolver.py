
import numpy as np
import sys, re

from scipy.optimize import root, minimize, newton


class FullSystemSolver(object):
    """FullSystemSolver provides a solver to determine minimum LCOE."""

    def __init__(self, system):
        self.full_system = system

    def minimizeFunctionBrownfield(self, initial_m_dot):

        initial_m_dot = initial_m_dot[0]

        if len(self.m_dots) > 10 and initial_m_dot >= (self.max_m_dot - 1e-8):
            print('NaN')
            return np.nan

        if len(self.m_dots) > 6 and np.isclose([np.mean(self.m_dots)], [self.initial_m_dot], rtol = 1e-2):
            print('STOP')
            return False
        else:
            self.m_dots.append(initial_m_dot)

        print('Trying a mass flowrate of %s'%initial_m_dot)

        try:
            system = self.full_system.solve(initial_m_dot, self.time_years)
            output = self.full_system.gatherOutput()
            output_val = output.capital_cost_model.LCOE_brownfield.LCOE

            if not np.isnan(output_val):
                self.test = np.array([initial_m_dot, output_val])

        except Exception as error:
            print(str(error))
            # regex = re.compile('Saturation pressure \[([+-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)) Pa\] corresponding to T \[([+-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)) K\] is within 1e-4 % of given p')
            # string = """Saturation pressure [7.2315e+06 Pa] corresponding to T [303.258 K] is within 1e-4 % of given p [7.2315e+06 Pa] : PropsSI("HMASS","P",7231504.074,"T",303.257922,"CO2")"""
            # re.match(regex, string)
            if str(error).find(':ExceedsMaxProductionPumpPressure') > -1:
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

        old_var = np.nan
        best_var = np.nan
        best_m_dot = np.nan

        #Starting Massflow Guess (add 10)
        m_dot_IP = 10

        peaks = 0
        d_m_dot = 10
        d_m_dot_multiplier = 0.21 #0.25

        point_towards_zero = False
        change_m_dot_sign = False

        # reversingFraction is the change past minimum to reverse
        reversing_fraction = 3e-2

        while peaks < 5:
            # Make sure mass flowrate is positive
            if m_dot_IP <= 0:
                d_m_dot = -1 * d_m_dot * d_m_dot_multiplier
                peaks = peaks + 1
                m_dot_IP = abs(d_m_dot)
                print('mdot<=0; peak %s'%peaks)

            print('Trying a mass flowrate of %s'%m_dot_IP)
            try:
                system = self.full_system.solve(m_dot_IP, time_years)
                output = self.full_system.gatherOutput()
                var_out = output.capital_cost_model.LCOE_brownfield.LCOE

                # See if this mass flowrate is better than previous values

                if np.isnan(var_out):
                    # make sure it is trending towards mdot of zero
                    point_towards_zero = True
                else:
                    if var_out < best_var or np.isnan(best_var):
                        best_var = var_out
                        best_m_dot = m_dot_IP

                    if peaks == 0:
                        # before first peak, exceed best val by
                        # reversingfraction
                        reversing_threshold = best_var * (1 + reversing_fraction)
                        if var_out > reversing_threshold:
                            change_m_dot_sign = True
                            print('Minimize first crossing; peak %s'%peaks)
                    else:
                        if var_out > old_var:
                            change_m_dot_sign = True
                            print('Minimize crossing; peak %s'%peaks)

                # set residual
                if np.isnan(var_out) or np.isnan(old_var) or d_m_dot == 0:
                    residual = 0
                else:
                    residual = abs(var_out - old_var) / abs(d_m_dot)

                # set current values as old value
                old_var = var_out

            except Exception as ex:
                raise ex
                print(str(ex))
                old_var = np.nan
                residual = 0
                point_towards_zero = True

            # # TODO: add raise exception

            if point_towards_zero:
                point_towards_zero = False
                if d_m_dot > 0:
                    change_m_dot_sign = True
                    print('Pointing to zero.')

            if change_m_dot_sign:
                change_m_dot_sign = False
                d_m_dot = -1 * d_m_dot * d_m_dot_multiplier
                peaks += 1
                print('Crossed peak %s, Changing direction.'%peaks)

            m_dot_IP = m_dot_IP + d_m_dot


        system = self.full_system.solve(m_dot_IP, time_years)
        output = self.full_system.gatherOutput()
        LCOE_brownfield = output.capital_cost_model.LCOE_brownfield.LCOE

        print(m_dot_IP, LCOE_brownfield)
        return m_dot_IP

    def minimizeLCOEBrownfieldAuto(self, time_years):

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
