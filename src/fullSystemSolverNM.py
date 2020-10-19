import numpy as np
import sys, re
from scipy.optimize import root, minimize, newton

from src.fullSystemSolverBase import FullSystemSolverBase

class FullSystemSolverNM(FullSystemSolverBase):
    """
    FullSystemSolverBase provides a solver to determine the optimum flow rate
    for a minimum or maximum of a given output variable.
    This solver uses the scipy minimize algorithm "Nelder-Mead".
    !!!!!!!!!!! This solver is not robust to local minima !!!!!!!!!!!
    """

    def __init__(self, system):
        super().__init__(system)

    def solve(self, time_years):

        print('You use the "Nelder-Mead" to find the optimum mdot \n'
              '!!!!!!!!!!! This solver is not robust to local minima !!!!!!!!!!!')

        self.initial_m_dot = 20.
        self.m_dots = []
        self.max_m_dot = 1e4
        self.time_years = time_years

        while True:

            sol = minimize(self.getOutputForMdot, (self.initial_m_dot),
            method='Nelder-Mead', tol=1e-3)

            if not sol.fun:
                self.initial_m_dot = self.initial_m_dot * 0.1
                self.m_dots = []
                if self.initial_m_dot < 1e-3:
                    return False
            else:
                return sol.x[0]
