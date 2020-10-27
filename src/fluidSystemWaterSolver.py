import numpy as np

from scipy.optimize import root, minimize, newton, brentq

from utils.fluidStateFromPT import FluidStateFromPT
from utils.solver import Solver


class FluidSystemWaterSolver(object):
    """FluidSystemWaterSolver provides a solver to determine water injection temperature."""

    def __init__(self, system):
        self.fluid_system = system

    def minimizeFunction(self, initialT, m_dot, time_years):
        # print('\n### new T ###')
        # print(initialT)
        # print(self.initial_P, initialT)
        initial_state = FluidStateFromPT(self.initial_P, initialT, self.fluid_system.fluid)
        system_state = self.fluid_system.solve(initial_state, m_dot, time_years)
        # self.initial_P =  self.fluid_system.pump.P_inj_surface
        # print(self.fluid_system.pump.P_inj_surface, system_state.T_C() - initialT, initialT)
        diff = (system_state.T_C() - initialT)
        # print(system_state.T_C())
        # print(diff)
        # print('### T done ###\n')
        print('find T ', system_state.T_C())
        return diff

    def minimizeFunctionOpt(self, initialT, m_dot, time_years):

        dT_inj = initialT
        solv = Solver()
        while abs(dT_inj) > 0.5:# 1e-3:

            initial_state = FluidStateFromPT(self.initial_P, initialT, self.fluid_system.fluid)
            system_state = self.fluid_system.solve(initial_state, m_dot, time_years)
            self.initial_P =  self.fluid_system.pump.P_inj_surface
            dT_inj = initialT - system_state.T_C()

            initialT = solv.addDataAndEstimate(initialT, dT_inj)

            if np.isnan(initialT):
                initialT = system_state.T_C()

            # add bounds
            if initialT < 1:
                initialT = 1

            # add upper bounds that make sense
            T_prod_surface_C = self.fluid_system.pump.well.results.T_C_f[-1]
            if initialT > T_prod_surface_C and initialT > 50:
                initialT = T_prod_surface_C


    def solve(self, m_dot, time_years):

        P_system_min = FluidStateFromPT.getPFromTQ(self.fluid_system.params.T_reservoir(), 0, self.fluid_system.fluid) + 1e5
        self.initial_P = P_system_min + 1e6

        state_T = 60.

        self.minimizeFunctionOpt(state_T, m_dot, time_years)

        # sol = newton(self.minimizeFunction, state_T, args=(m_dot, time_years),
        # maxiter = 50, tol = 1.e-14, rtol = 1e-8)
        # print('final T ', sol)
        # 1/0

        # try:
        #     sol = newton(self.minimizeFunction, state_T, args=(m_dot, time_years),
        #     maxiter = 50, tol = 1.e-14, rtol = 1e-8)
        #     # print('P_Pa, ', self.initial_P)
        #     # print('sol,  ', sol)
        # except Exception as error:
        #     if str(error).find('SemiAnalyticalWell:BelowSaturationPressure') > -1:
        #         print(str(error))
        #         sol = newton(self.minimizeFunction, 40., args=(m_dot, time_years),
        #         maxiter = 50, tol = 1.e-14, rtol = 1e-8)
        #     else:
        #         raise(error)


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

            # sol = brentq(self.minimizeFunction, 40., 1e-6, args=(m_dot, time_years),
            # maxiter = 25, full_output = True, xtol = 1.e-4, rtol = 1e-8)

            # sol = root(self.minimizeFunction, state_T, args=(m_dot, time_years),
            # method ='lm', tol = 1.e-4, options = {'xtol' : 1e-3})

        # initial_state = FluidStateFromPT(self.initial_P, sol, self.fluid_system.fluid)
        # system_state = self.fluid_system.solve(initial_state, m_dot, time_years)
        # diff = abs(system_state.T_C() - sol)
        # if diff > 1e-5:
        #     print('FluidSystemWater:NotConverged - '
        #     'Difference between injection and power plant output T is ', diff)


    def gatherOutput(self):
        return self.fluid_system.gatherOutput()
