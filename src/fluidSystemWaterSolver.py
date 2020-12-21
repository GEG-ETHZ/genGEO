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

from scipy.optimize import root, minimize, newton, brentq

from utils.fluidState import FluidState
from utils.solver import Solver


class FluidSystemWaterSolver(object):
    """FluidSystemWaterSolver provides a solver to determine water injection temperature."""

    def __init__(self, system):
        self.fluid_system = system

    def solve(self):

        self.fluid_system.params.initial_P = 1e6 + self.fluid_system.params.P_system_min()
        self.fluid_system.params.initial_T = 60.
        initial_T = self.fluid_system.params.initial_T

        dT_inj = np.nan
        dT_loops = 1
        solv = Solver()
        while np.isnan(dT_inj) or abs(dT_inj) >= 0.5:

            initial_state = FluidState.getStateFromPT(self.fluid_system.params.initial_P, initial_T, self.fluid_system.params.working_fluid)
            system_state = self.fluid_system.solve(initial_state)

            dT_inj = initial_T - system_state.pp.state.T_C

            initial_T = solv.addDataAndEstimate(initial_T, dT_inj)

            if np.isnan(initial_T):
                initial_T = system_state.pp.state.T_C

            # add lower bounds
            if initial_T < 1:
                initial_T = 1

            # add upper bounds
            T_prod_surface_C = system_state.pump.well.state.T_C
            if initial_T > T_prod_surface_C and initial_T > 50:
                initial_T = T_prod_surface_C

            if dT_loops >  10:
                print('GenGeo::Warning:FluidSystemWaterSolver:dT_loops is large: %s'%dT_loops)
            dT_loops += 1

        # check if silica precipitation is allowed
        if self.fluid_system.params.silica_precipitation:
            # prevent silica precipitation by DiPippo 1985
            maxSurface_dT = 89
            if (T_prod_surface_C - initial_T) > maxSurface_dT:
                raise Exception('GenGeo::FluidSystemWaterSolver:ExceedsMaxTemperatureDecrease - '
                            'Exceeds Max Temp Decrease of  %.3f C to prevent silica precipitation!'%(maxSurface_dT))
        else:
            maxSurface_dT = np.inf

        return system_state

    def minimizeFunction(self, initial_T):

        initial_state = FluidState.getStateFromPT(self.initial_P, initial_T, self.fluid_system.fluid)
        system_state = self.fluid_system.solve(initial_state)

        diff = (system_state.pp.state.T_C - initial_T)

        print('find T ', system_state.pp.state.T_C)
        return diff

    def solveMinimize(self):

        self.initial_P = 1e6 +  self.fluid_system.params.P_system_min()
        initial_T = 60.

        try:
            sol = newton(self.minimizeFunction, initial_T,
            maxiter = 50, tol = 1.e-14, rtol = 1e-8)

        except Exception as ex:
            if str(ex).find('BelowSaturationPressure') > -1:
                print(str(ex))
                sol = newton(self.minimizeFunction, 40.,
                maxiter = 50, tol = 1.e-14, rtol = 1e-8)
            else:
                raise ex


        # state_T = 1.

        # while state_T < 100.:
        #     try:
        #         sol = newton(self.minimizeFunction, state_T, args=(m_dot, time_years), #full_output = True,
        #         maxiter = 25, tol = 1.e-4, rtol = 1e-8)
        #         state_T = 111.
        #     except Exception as error:
        #         state_T = state_T + 20.
        #         if state_T >= 100.:
        #             raise(error)

        # state_T = 100.
        # while state_T > 0.:
        #     try:
        #         sol = newton(self.minimizeFunction, state_T, args=(m_dot, time_years),
        #         maxiter = 25, full_output = True, tol = 1.e-4, rtol = 1e-8)
        #         state_T = -1.
        #     except Exception as error:
        #         state_T = state_T - 20.
        #         if state_T <= 0.:
        #             raise(error)


        # initial_state = FluidStateFromPT(self.initial_P, sol, self.fluid_system.fluid)
        # system_state = self.fluid_system.solve(initial_state, m_dot, time_years)
        # diff = abs(system_state.T_C() - sol)
        # if diff > 1e-5:
        #     print('FluidSystemWater:NotConverged - '
        #     'Difference between injection and power plant output T is ', diff)
