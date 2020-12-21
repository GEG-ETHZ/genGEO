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

from src.semiAnalyticalWellResults import SemiAnalyticalWellResults

from utils.constantsAndPaths import ConversionConstants
from utils.fluidState import FluidState
from utils.frictionFactor import frictionFactor
from models.simulationParameters import SimulationParameters

class SemiAnalyticalWell(object):
    """SemiAnalyticalWell to compute heat transport with fluid flow in a well
        and analytical conduction in the surrounding rock."""

    def __init__(self, params = None, T_e_initial = 15., dz_total = 0., dr_total = 0., m_dot_multiplier = 1, **kwargs):
        self.params = params
        if self.params == None:
            self.params = SimulationParameters(**kwargs)
        # dimensions
        # dz_total is the change in elevation, negative if going into the earth
        self.dz_total = dz_total
        # dr_total is the deviation from vertical, always positive
        self.dr_total = dr_total
        # T_e_initial is the beginning rock temp
        # default is set to be 15 C
        self.T_e_initial = T_e_initial
        # m_dot_multiplier is the multiplier applied to the mass flowrate within the wells
        # default is 1
        self.m_dot_multiplier = m_dot_multiplier

    def solve(self, initial_state):

        m_dot = self.params.m_dot_IP * self.m_dot_multiplier

        # results
        results = SemiAnalyticalWellResults(self.params.well_segments, self.params.working_fluid)

        P_f_initial = initial_state.P_Pa
        T_f_initial = initial_state.T_C
        time_seconds = self.params.time_years * ConversionConstants.secPerYear

        # set geometry
        dz = self.dz_total/self.params.well_segments             # m
        dr = self.dr_total/self.params.well_segments             # m
        dL = (dz**2 + dr**2)**0.5                       # m
        A_c = np.pi * self.params.well_radius**2        # m**2
        P_c = np.pi * 2 * self.params.well_radius       # m

        # set states
        results.T_C_f[0] = T_f_initial             # C
        results.T_C_e[0] = self.T_e_initial        # C
        results.P_Pa[0] = P_f_initial              # Pa
        results.h_Jkg[0] = FluidState.getStateFromPT(results.P_Pa[0], results.T_C_f[0], self.params.working_fluid).h_Jkg
        results.rho_kgm3[0] = FluidState.getStateFromPT(results.P_Pa[0], results.T_C_f[0], self.params.working_fluid).rho_kgm3
        results.v_ms[0] = m_dot / A_c / results.rho_kgm3[0]  #m/s

        # Calculate the Friction Factor
        # Use Colebrook-white equation for wellbore friction loss.
        # Calculate for first element and assume constant in remainder
        ff = frictionFactor(self.params.well_radius, results.P_Pa[0], results.h_Jkg[0], \
                            m_dot, self.params.working_fluid, self.params.epsilon)

        alpha_rock = self.params.k_rock/self.params.rho_rock/self.params.c_rock  #D rock
        t_d = alpha_rock*time_seconds/(self.params.well_radius**2)  #dim
        if t_d < 2.8:
            beta = ((np.pi*t_d)**-0.5 + 0.5 - 0.25*(t_d/np.pi)**0.5 + 0.125*t_d)
        else:
            beta = (2/(np.log(4*t_d)-2*0.58) - 2*0.58/(np.log(4*t_d)-2*0.58)**2)

        # loop over all well segments
        for i in range(1, self.params.well_segments+1):
            results.z_m[i] = results.z_m[i-1] + dz

            # far-field rock temp
            results.T_C_e[i] = results.T_C_e[0] - results.z_m[i] * self.params.dT_dz
            # fluid velocity
            results.v_ms[i] = m_dot / A_c / results.rho_kgm3[i-1]  #m/s

            # Calculate Pressure
            results.delta_P_loss[i] = ff * dL / ( 2 * self.params.well_radius) * \
                                            results.rho_kgm3[i-1] * \
                                            results.v_ms[i]**2. / 2.  #Pa
            results.rho_kgm3[i] = FluidState.getStateFromPh(results.P_Pa[i-1], results.h_Jkg[i-1], self.params.working_fluid).rho_kgm3
            results.P_Pa[i]     = results.P_Pa[i-1] - results.rho_kgm3[i] * \
                                        self.params.g * dz - results.delta_P_loss[i]

            # Throw exception if below saturation pressure of water at previous temperature
            if self.params.working_fluid.lower() == 'water':
                P_sat = FluidState.getStateFromTQ(results.T_C_f[i-1], 0, self.params.working_fluid).P_Pa
                if results.P_Pa[i] < P_sat:
                    raise Exception('GenGeo::SemiAnalyticalWell:BelowSaturationPressure - '
                    'Below saturation pressure of water at %s m !' %(results.z_m[i]))

            h_noHX = results.h_Jkg[i-1] - self.params.g * dz
            T_noHX = FluidState.getStateFromPh(results.P_Pa[i], h_noHX, self.params.working_fluid).T_C
            results.cp_JK[i] = FluidState.getStateFromPh(results.P_Pa[i], h_noHX, self.params.working_fluid).cp_JK

            #Find Fluid Temp
            if not self.params.useWellboreHeatLoss:
                results.T_C_f[i] = T_noHX
                results.h_Jkg[i] = h_noHX
                results.q[i] = 0.
            else:
                # See Zhang, Pan, Pruess, Finsterle (2011). A time-convolution
                # approach for modeling heat exchange between a wellbore and
                # surrounding formation. Geothermics 40, 261-266.
                x = dL * P_c * self.params.k_rock * beta / self.params.well_radius
                y = m_dot * results.cp_JK[i]
                if math.isinf(x):
                    results.T_C_f[i] = results.T_C_e[i]
                else:
                    results.T_C_f[i] = (y * T_noHX + x *results. T_C_e[i]) / (x + y)
                results.q[i] = y * (T_noHX - results.T_C_f[i])
                results.h_Jkg[i] = FluidState.getStateFromPT(results.P_Pa[i], results.T_C_f[i], self.params.working_fluid).h_Jkg
        # make sure state object is set
        results.createFinalState()

        return results
