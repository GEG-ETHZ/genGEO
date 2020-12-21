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

from  utils.fluidState import FluidState

def frictionFactor(well_radius, P, h, m_dot, fluid, epsilon):

    if well_radius == None or P == None or h == None or m_dot == None or fluid == None or epsilon == None:
        return 0

    rho_fluid = FluidState.getStateFromPh(P, h, fluid).rho_kgm3
    mu = FluidState.getStateFromPh(P, h, fluid).mu_Pas

    A_c = np.pi * well_radius**2
    V = m_dot / A_c / rho_fluid

    # Relative Roughness
    rr = epsilon / (2 * well_radius)

    # Reynolds Number
    Re = rho_fluid * V * (2*well_radius) / mu

    # Use Haaland (1983) Approximation
    c1 = -1.8 * np.log10((6.9 / Re) + (rr / 3.7)**1.11)
    return (1 / c1)**2
