# Licensed under LGPL 2.1, please see LICENSE for details
# https://www.gnu.org/licenses/lgpl-2.1.html
#
# The work on this project has been performed at the GEG Group at ETH Zurich:
# --> https://geg.ethz.ch
#
# The initial version of this file has been implemented by:
#
#     Philipp Schaedle (https://github.com/philippschaedle)
#     Benjamin M. Adams
#
# Further changes are done by:
#

############################
import numpy as np
import re

from src.fullSystemOptMdotBase import FullSystemOptMdotBase

class FullSystemSolverBase(FullSystemOptMdotBase):
    """
    FullSystemSolverBase provides a solver to solve
    for a given flow rate.
    """
    def __init__(self, system):
        super().__init__(system)
        self.m_dots = []
        self.max_m_dot = np.inf

    def getOutputForMdot(self, initial_m_dot):

        initial_m_dot = initial_m_dot[0]

        if len(self.m_dots) > 10 and initial_m_dot >= (self.max_m_dot - 1e-8):
            return np.nan

        if len(self.m_dots) > 6 and np.isclose([np.mean(self.m_dots)], [self.initial_m_dot], rtol = 1e-2):
            return False
        else:
            self.m_dots.append(initial_m_dot)

        print('Trying a mass flowrate of %s'%initial_m_dot)

        try:
            system = self.full_system.solve(initial_m_dot, self.time_years)
            output_val = self.getTargetVar(system)

            # if not np.isnan(output_val):
            #     self.test = np.array([initial_m_dot, output_val])

        except Exception as ex:
            print(str(ex))
            # regex = re.compile('Saturation pressure \[([+-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)) Pa\] corresponding to T \[([+-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)) K\] is within 1e-4 % of given p')
            # string = """Saturation pressure [7.2315e+06 Pa] corresponding to T [303.258 K] is within 1e-4 % of given p [7.2315e+06 Pa] : PropsSI("HMASS","P",7231504.074,"T",303.257922,"CO2")"""
            # re.match(regex, string)
            if str(ex).find(':ExceedsMaxProductionPumpPressure') > -1 \
            or str(ex).find(':BelowSaturationPressure') > -1:
                self.max_m_dot = initial_m_dot
                output_val = np.nan
            # elif str(ex).find('FluidSystemCO2:TurbinePowerNegative') > -1 \
            # or str(ex).find('Saturation pressure ') > -1: # # TODO: fix this to regular expression ex triggered in semiAnalyticalWell.py line: 123
            #     n = 2
            #     while n <= 3:
            #         try:
            #             diff = initial_m_dot[0] - self.test[0]
            #             tmp_m_dot = self.test[0] + n * diff
            #             system = self.full_system.solve(tmp_m_dot, self.time_years)
            #             output_val_tmp = self.getTargetVar(system)
            #             output_val = np.interp(initial_m_dot, [self.test[0], tmp_m_dot], [self.test[1], output_val_tmp])[0]
            #             n = 11
            #         except:
            #             n += 1
            #             output_val = np.nan
            else:
                raise ex

        return output_val
