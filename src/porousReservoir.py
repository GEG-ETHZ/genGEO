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
import math
import numpy as np

from src.porousReservoirBase import PorousReservoirBase
from src.porousReservoirResults import PorousReservoirResults

from utils.fluidState import FluidState
from utils.depletionCurve import depletionCurve
from utils.constantsAndPaths import ConversionConstants

from models.simulationParameters import SimulationParameters
from models.wellFieldType import WellFieldType

class PorousReservoir(PorousReservoirBase):
    """PorousReservoir holds methods to compute reservoir heat extraction.
    """
    def __init__(self, params = None,
            modelPressureTransient = False,
            modelTemperatureDepletion = True, **kwargs):
        self.params = params
        if self.params == None:
            self.params = SimulationParameters(**kwargs)
        self.modelPressureTransient = modelPressureTransient
        self.modelTemperatureDepletion = modelTemperatureDepletion
        super().__init__()

    def solve(self, initialState):
        results = PorousReservoirResults()

        A_reservoir = (self.params.well_spacing**2)/2.

        # # TODO: fix this. using 365 instead of 365.25 changes the results
        time_seconds = self.params.time_years * ConversionConstants.secPerYear

        # from output properties
        prod_State = FluidState.getStateFromPT(self.params.P_reservoir(), self.params.T_reservoir(), self.params.working_fluid)

        mu_fluid = (initialState.mu_Pas + prod_State.mu_Pas)/2
        rho_fluid = (initialState.rho_kgm3 + prod_State.rho_kgm3)/2
        cp_fluid = (initialState.cp_JK + prod_State.cp_JK)/2


        # reservoir impedance
        if self.params.well_spacing <= self.params.well_radius:
            RI = 0
        else:
            if self.params.wellFieldType == WellFieldType._5Spot_SharedNeighbor \
                    or self.params.wellFieldType == WellFieldType._5Spot \
                    or self.params.wellFieldType == WellFieldType._5Spot_Many:
                A_c_rock = np.log( (4 * self.params.well_spacing) / (2 * self.params.well_radius * np.pi ))
                RI = mu_fluid/rho_fluid/self.params.transmissivity * A_c_rock
            elif self.params.wellFieldType == WellFieldType.Doublet:
                A_c_rock = np.log( (self.params.well_spacing) / (2*self.params.well_radius) )
                RI = mu_fluid/rho_fluid/self.params.transmissivity/np.pi * A_c_rock
            else:
                raise Exception('GenGeo::PorousReservoir:unknownReservoirConfiguration - '
                'Unknown Reservoir Configuration')

        # Calculate heat extracted
        dT_initial = self.params.T_reservoir() - initialState.T_C
        res_energy = A_reservoir * self.params.reservoir_thickness * self.params.rho_rock * self.params.c_rock * dT_initial

        # Model pressure transient (Figure 4.2, Adams (2015)), only for CO2 drying out
        if self.modelPressureTransient == True and self.params.working_fluid.lower()=='co2':
            R = self.params.well_spacing
            nu_inj_fluid = initialState.mu_Pas / initialState.rho_kgm3

            tau = time_seconds * nu_inj_fluid / R**2.

            b = -0.4574 * (abs(R)/abs(self.params.depth))**-0.27
            b = np.clip(b, -0.9, -0.4)
            a = 1.0342 * np.exp(10.989*b)
            a = np.clip(a, 0., 0.012)
            coeff = (a*tau**b + 1)
            RI = RI * coeff


        # Model temperature transient (Figure 4.5, Adams (2015))
        # First-principles model for all fluids
        # At time zero, output is the same as no temp depletion
        if self.modelTemperatureDepletion == True and self.params.time_years > 0.:
            # p1
            coeff1 = self.params.reservoir_thickness * self.params.k_rock * initialState.rho_kgm3 / self.params.m_dot_IP / self.params.rho_rock / self.params.c_rock
            p1 = -890.97*coeff1 - 1.30
            # limit to regression
            p1 = np.minimum(p1, -1.3)
            # p2
            p2 = 0.4095*np.exp(-0.7815 * p1)
            # p3
            R = self.params.well_spacing
            coeff2 = self.params.k_rock * R * initialState.rho_kgm3 / self.params.rho_rock / initialState.cp_JK / self.params.m_dot_IP
            # A more realistic, exponential relation is found by fitting the same data
            # with an exponential, instead of linear, curve.
            p3 = 1 - (1.4646 * np.exp(-377.3*coeff2))
            # p3 cant be greater than 1 or less than 0.
            p3 = np.clip(p3, 0., 1.)
            # psi
            # numerically integrate to get heat depletion
            # have roughly one increment a year
            increments = np.maximum(np.round(self.params.time_years), 1.)
            dt = time_seconds / increments
            res_energy_extracted = 0
            gamma = 1
            for i in range(int(increments)):
                res_energy_extracted = self.params.m_dot_IP * cp_fluid * gamma * dT_initial * dt + res_energy_extracted
                results.psi = res_energy_extracted / res_energy
                gamma = depletionCurve(results.psi, p1, p2, p3)
                # last gamma is res temp
        else:
            # gamma of 1 is undepleted reservoir
            gamma = 1
            res_energy_extracted = self.params.m_dot_IP * cp_fluid * gamma * dT_initial * time_seconds
            if res_energy == 0:
                results.psi = 0
            else:
                results.psi = res_energy_extracted / res_energy


        results.dP = self.params.m_dot_IP * RI
        end_P_Pa = initialState.P_Pa - results.dP
        end_T_C = gamma * dT_initial + initialState.T_C
        results.state = FluidState.getStateFromPT(end_P_Pa, end_T_C, self.params.working_fluid)
        results.heat = self.params.m_dot_IP * (results.state.h_Jkg - initialState.h_Jkg)

        return results
