
from scipy.optimize import root, minimize, newton

from utils.fluidStateFromPT import FluidStateFromPT


class FluidSystemWaterSolver(object):
    """FluidSystemWaterSolver provides a solver to determine water injection temperature."""

    def __init__(self, system):
        self.fluid_system = system

    def minimizeFunction(self, initialT):
        initial_state = FluidStateFromPT(1.e6, initialT, self.fluid_system.fluid)
        system_state = self.fluid_system.solve(initial_state, self.m_dot, self.time_years)
        return (system_state.T_C() - initialT)

    def solve(self, m_dot, time_years):

        self.m_dot = m_dot
        self.time_years = time_years

        state_T = 60.
        # optionsSolve={}
        # optionsSolve['xtol'] = 1e-3
        # sol = root(self.minimizeFunction, state_T, method ='df-sane')#, options = optionsSolve)
        # sol = minimize(self.minimizeFunction, (state_T), method='dogleg', tol=1e-20)

        sol = newton(self.minimizeFunction, state_T, tol =  1., rtol = 1e-3)
