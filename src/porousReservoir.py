import math
import numpy as np

from utils.fluidStates import FluidState
from utils.depletionCurve import depletionCurve
from utils.globalConstants import globalConstants
from src.porousReservoirResults import PorousReservoirResults

class PorousReservoir(object):
    """PorousReservoir holds methods to compute reservoir heat extraction.
    """
    def __init__(self, params, well_spacing, permeability,
            thickness, depth, dT_dz, wellRadius, T_surface_rock, fluid,
            reservoirConfiguration,
            modelPressureTransient,
            modelTemperatureDepletion):
        self.params = params
        self.well_spacing = well_spacing
        self.transmissivity = permeability * thickness
        self.thickness =  thickness
        self.T_surface_rock = T_surface_rock
        self.depth = depth
        self.dT_dz = dT_dz
        self.wellRadius = wellRadius
        self.fluid = fluid
        self.reservoirT = abs(self.depth) * abs(self.dT_dz) + self.T_surface_rock
        self.P_reservoir = self.depth * 1000. * self.params.g
        self.reservoirConfiguration = reservoirConfiguration
        self.modelPressureTransient = modelPressureTransient
        self.modelTemperatureDepletion = modelTemperatureDepletion

    def solve(self, initialState, m_dot, time_years):
        results = PorousReservoirResults(self.fluid, self.well_spacing)

        self.m_dot =  m_dot
        self.time_years = time_years
        self.initialState = initialState

        A_reservoir = (self.well_spacing**2)/2.

        # # TODO: fix this. using 365 instead of 365.25 changes the results
        time_seconds = 3600 * 24 * 365 * self.time_years #* globalConstants.secPerYear

        self.depth = abs(self.depth)

        # from output properties
        prod_State = FluidState.getStateFromPT(self.P_reservoir, self.reservoirT, self.fluid)

        mu_fluid = (self.initialState.v_Pas + prod_State.v_Pas)/2
        rho_fluid = (self.initialState.rho_kgm3 + prod_State.rho_kgm3)/2
        cp_fluid = (self.initialState.cp_JK + prod_State.cp_JK)/2


        # reservoir impedance
        if self.well_spacing <= self.wellRadius:
            self.RI = 0
        else:
            if self.reservoirConfiguration == '5spot':
                A_c_rock = np.log( (4 * self.well_spacing) / (2 * self.wellRadius * np.pi ))
                self.RI = mu_fluid/rho_fluid/self.transmissivity * A_c_rock
            elif self.reservoirConfiguration == 'Doublet':
                A_c_rock = np.log( (self.well_spacing) / (2*self.wellRadius) )
                self.RI = mu_fluid/rho_fluid/self.transmissivity/np.pi * A_c_rock
            else:
                raise ValueError('PorousReservoir:unknownReservoirConfiguration',
                'Unknown Reservoir Configuration')

        # Calculate heat extracted
        dT_initial = self.reservoirT - self.initialState.T_C
        res_energy = A_reservoir * self.thickness * self.params.rho_rock * self.params.c_rock * dT_initial

        # Model pressure transient (Figure 4.2, Adams (2015)), only for CO2 drying out
        if self.modelPressureTransient == True and self.fluid=='CO2':
            R = self.well_spacing
            nu_inj_fluid = self.initialState.v_Pas / self.initialState.rho_kgm3
            # fit doesn't work before 2 years
            if self.time_years < 2.:
                tau = (globalConstants.secPerYear * 2.) * nu_inj_fluid / R**2.
            else:
                tau = time_seconds * nu_inj_fluid / R**2.

            b = -0.4574 * (abs(R)/abs(self.depth))**-0.27
            b = np.clip(b, -0.9, -0.4)
            a = 1.0342 * np.exp(10.989*b)
            a = np.clip(a, 0., 0.012)
            coeff = (a*tau**b + 1)
            self.RI = self.RI * coeff


        # Model temperature transient (Figure 4.5, Adams (2015))
        # First-principles model for all fluids
        # At time zero, output is the same as no temp depletion
        if self.modelTemperatureDepletion == True and self.time_years > 0.:
            # p1
            coeff1 = self.thickness * self.params.k_rock * self.initialState.rho_kgm3 / self.m_dot / self.params.rho_rock / self.params.c_rock
            p1 = -890.97*coeff1 - 1.30
            # limit to regression
            p1 = np.minimum(p1, -1.3)
            # p2
            p2 = 0.4095*np.exp(-0.7815 * p1)
            # p3
            R = self.well_spacing
            coeff2 = self.params.k_rock * R * self.initialState.rho_kgm3 / self.params.rho_rock / self.initialState.cp_JK / self.m_dot
            # A more realistic, exponential relation is found by fitting the same data
            # with an exponential, instead of linear, curve.
            p3 = 1 - (1.4646 * np.exp(-377.3*coeff2))
            # p3 cant be greater than 1 or less than 0.
            p3 = np.clip(p3, 0., 1.)
            # psi
            # numerically integrate to get heat depletion
            # have roughly one increment a year
            increments = np.maximum(np.round(self.time_years), 1.)
            dt = time_seconds / increments
            res_energy_extracted = 0
            gamma = 1
            for i in range(int(increments)):
                res_energy_extracted = self.m_dot * cp_fluid * gamma * dT_initial * dt + res_energy_extracted
                self.psi = res_energy_extracted / res_energy
                gamma = depletionCurve(self.psi, p1, p2, p3)
                # last gamma is res temp
        else:
            # gamma of 1 is undepleted reservoir
            gamma = 1
            res_energy_extracted = self.m_dot * cp_fluid * gamma * dT_initial * time_seconds
            if res_energy == 0:
                self.psi = 0
            else:
                self.psi = res_energy_extracted / res_energy

        self.end_P_Pa = self.initialState.P_Pa - (self.m_dot * self.RI)
        self.end_T_C = gamma * dT_initial + self.initialState.T_C
        return FluidState.getStateFromPT(self.end_P_Pa, self.end_T_C, self.fluid)

    def gatherOutput(self):
        output = PorousReservoirResults(self.fluid, self.well_spacing)
        # Calculate final values
        output.dP = self.m_dot * self.RI
        output.end_P_Pa = self.initialState.P_Pa - output.dP
        output.end_T_C = self.end_T_C
        output.end_h_Jkg = FluidState.getHFromPT(self.end_P_Pa, self.end_T_C, self.fluid)
        output.heat = self.m_dot * (output.end_h_Jkg - self.initialState.h_Jkg)
        output.psi =  self.psi
        output.reservoirT = self.reservoirT
        return output
