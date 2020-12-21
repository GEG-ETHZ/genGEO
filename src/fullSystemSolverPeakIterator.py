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

from src.fullSystemOptMdotBase import FullSystemOptMdotBase

from enum import Enum

class FullSystemSolverPeakIterator(FullSystemOptMdotBase):
    """
    FullSystemSolverPeakIterator implements a ...
    """

    def __init__(self, system):
        super().__init__(system)

    def getDirection(self):
        return SolverOptimizationType.Maximize

    def getTargetVar(self):
        return None

    def solve(self):
        old_var = np.nan
        best_var = np.nan
        #best_m_dot = np.nan

        #Starting Massflow Guess (add 10)
        m_dot_IP = 10

        peaks = 0
        d_m_dot = 10
        d_m_dot_multiplier = 0.21 #0.25

        point_towards_zero = False
        change_m_dot_sign = False

        # reversingFraction is the change past minimum to reverse
        reversing_fraction = 3e-2

        opt_dir = np.sign(self.getDirection().value)
        opt_dir_dict = {-1: 'Minimize',
                         1: 'Maximize'}

        while peaks < 5:
            # Make sure mass flowrate is positive
            if m_dot_IP <= 0:
                d_m_dot = -1 * d_m_dot * d_m_dot_multiplier
                peaks = peaks + 1
                m_dot_IP = abs(d_m_dot)
                print('mdot<=0; peak %s'%peaks)

            print('Trying a mass flowrate of %.2f' %m_dot_IP)
            try:
                self.full_system.params.m_dot_IP = m_dot_IP
                solveResult = self.full_system.solve()
                var_out = self.getTargetVar(solveResult)

                # See if this mass flowrate is better than previous values

                if np.isnan(var_out):
                    # make sure it is trending towards mdot of zero
                    point_towards_zero = True
                else:
                    # if var_out < best_var or np.isnan(best_var):
                    if np.sign(var_out - best_var) == opt_dir or np.isnan(best_var):
                        best_var = var_out
                        best_m_dot = m_dot_IP

                    if peaks == 0:
                        # before first peak, exceed best val by
                        # reversing fraction
                        reversing_threshold = best_var * (1 - reversing_fraction * opt_dir)
                        # if var_out > reversing_threshold:
                        if np.sign(reversing_threshold - var_out) == opt_dir:
                            change_m_dot_sign = True
                            print('%s first crossing; peak %s'%(opt_dir_dict[opt_dir],peaks))
                    else:
                        # if var_out > old_var:
                        if np.sign(old_var - var_out) == opt_dir:
                            change_m_dot_sign = True
                            print('%s crossing; peak %s'%(opt_dir_dict[opt_dir],peaks))

                # set residual
                if np.isnan(var_out) or np.isnan(old_var) or d_m_dot == 0:
                    residual = 0
                else:
                    residual = abs(var_out - old_var) / abs(d_m_dot)

                # set current values as old value
                old_var = var_out

            except Exception as ex:
                if str(ex).find(':ExceedsMaxProductionPumpPressure') > -1 \
                or str(ex).find(':BelowSaturationPressure') > -1:
                    print(str(ex))
                    old_var = np.nan
                    residual = 0
                    point_towards_zero = True
                else:
                    raise ex

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

        # one final solve with the optimum flowrate
        self.full_system.params.m_dot_IP = m_dot_IP
        results = self.full_system.solve()
        results.optMdot = m_dot_IP

        return results


class SolverOptimizationType(Enum):
    Minimize = -1
    Maximize = 1