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
