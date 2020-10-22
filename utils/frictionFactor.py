import numpy as np

from  utils.fluidStates import FluidState

def frictionFactor(well_radius, P, h, m_dot, fluid, epsilon):

    rho_fluid = FluidState.getRhoFromPh(P, h, fluid)
    mu = FluidState.getMuFromPh(P, h, fluid)

    A_c = np.pi * well_radius**2
    V = m_dot / A_c / rho_fluid

    # Relative Roughness
    rr = epsilon / (2 * well_radius)

    # Reynolds Number
    Re = rho_fluid * V * (2*well_radius) / mu

    # Use Haaland (1983) Approximation
    c1 = -1.8 * np.log10((6.9 / Re) + (rr / 3.7)**1.11)
    return (1 / c1)**2
