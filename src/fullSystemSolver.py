
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

        LCOE_old = np.nan
        LCOE = 0

        bestLCOE = 1e12
        bestLCOEMdot = np.nan

        #Starting Massflow Guess (add 10)
        m_dot_IP = 0

        peaks = 0
        d_m_dot = 10
        d_power_d_m_dot = 0
        d_m_dot_multiplier = 0.21 #0.25

        while peaks < 5:
            # if step is too big and crosses zero, make it smaller
            if m_dot_IP + d_m_dot <= 0:
                d_m_dot = -1 * d_m_dot * d_m_dot_multiplier
                peaks = peaks + 1
                print('mdot<=0; peak %s'%peaks)
                # get it to go back
                # W_net_IP_old = np.nan
                LCOE_old = np.nan
            m_dot_IP = m_dot_IP + d_m_dot
            print('Trying a mass flowrate of %s'%m_dot_IP)
            try:
                system = self.full_system.solve(m_dot_IP, time_years)
                output = self.full_system.gatherOutput()
                LCOE_brownfield = output.capital_cost_model.LCOE_brownfield.LCOE
                print('Done')
                # # should never be zero
                # if d_m_dot != 0 and not np.isnan(W_net_IP) and not np.isnan(W_net_IP_old):
                #     d_power_d_m_dot = abs((W_net_IP - W_net_IP_old)*1e6 / d_m_dot)

                print(LCOE_old)

                # See if LCOE has inflected
                # LCOE is NaN only at zero flowrate
                if np.isnan(LCOE_brownfield):
                    # If NaN LCOE, make sure it is trending towards mdot of
                    # zero
                    if d_m_dot > 0:
                        d_m_dot = -1 * d_m_dot * d_m_dot_multiplier
                        peaks = peaks + 1
                        print('LCOE_brownfield=NaN; peak %s'%peaks)

                elif not np.isnan(LCOE_old):
                    # only make change if both LCOE and LCOE_old are not
                    # nan
                    if LCOE_old < LCOE_brownfield:
                        d_m_dot = -1 * d_m_dot * d_m_dot_multiplier
                        peaks = peaks + 1
                        print('LCOE_old<LCOE_brownfield; peak %s'%peaks)

                LCOE = LCOE_brownfield


                # Save best values
                if not np.isnan(LCOE_brownfield) and LCOE_brownfield < bestLCOE:
                    bestLCOE = LCOE
                    bestLCOEMdot = m_dot_IP

                # Set current values to old values
                LCOE_old = LCOE

            except Exception as ex:
                print(str(ex))
                if d_m_dot > 0:
                    d_m_dot = -1 * d_m_dot * d_m_dot_multiplier
                    peaks = peaks + 1
                    print('exception; peak %s'%peaks)
                    LCOE_old = np.nan
                    LCOE_brownfield = np.nan
                    LCOE = np.nan
            # # TODO: add raise exception

        m_dot_IP = bestLCOEMdot
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
