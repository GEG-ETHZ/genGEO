import numpy as np

from scipy.optimize import root, minimize, newton

class FullSystemSolverBase(object):
    """FullSystemSolverBase provides a solver to determine the optimum flow rate for a minimum or maximum of a given output variable."""

    def __init__(self, system):
        self.full_system = system

    def getTargetVar(self):
        raise Exception('GenGeo::no target variable provided to find Minimum '
                        ' or Maximum of!')

    def getDirection(self):
        raise Exception('GenGeo::no direction provided to find Minimum '
                        ' or Maximum getDirection must be '
                        '-1 for minimum or'
                        ' 1 for maximum!')

    def findOptMdotForExtrema(self, time_years):

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

        opt_dir = np.sign(self.getDirection())
        opt_dir_dict = {-1: 'Minimize',
                         1: 'Maximize'}

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
                var_out = self.getTargetVar()

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
                print(str(ex))
                old_var = np.nan
                residual = 0
                point_towards_zero = True

            # # TODO: add raise exception
            # raise ex

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
        system = self.full_system.solve(m_dot_IP, time_years)

        t_val = self.getTargetVar()

        return m_dot_IP
