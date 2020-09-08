import math
import numpy as np

from utils.fluidStates import FluidState
from src.semiAnalyticalWellResults import SemiAnalyticalWellResults

class SemiAnalyticalWell(object):
    """SemiAnalyticalWell to compute heat transport with fluid flow in a well
        and analytical conduction in the surrounding rock."""

    def __init__(self, params, wellRadius, fluid, epsilon, dT_dz, dz_total = 0., dr_total = 0.):
        self.params = params
        # dimensions
        self.wellRadius = wellRadius
        self.dz_total = dz_total
        self.dr_total = dr_total
        # states
        self.fluid = fluid
        # model
        self.epsilon = epsilon
        self.dT_dz = dT_dz

    def solve(self, P_f_initial, T_f_initial, T_e_initial, time_seconds, m_dot):
        results = SemiAnalyticalWellResults(self.params.N_dx, self.fluid)

        # set geometry
        dz = self.dz_total/self.params.N_dx           # m
        dr = self.dr_total/self.params.N_dx           # m
        dL = (dz**2 + dr**2)**0.5    # m
        A_c = np.pi * self.wellRadius**2       # m**2
        P_c = np.pi * 2 * self.wellRadius      # m

        # set states
        results.T_C_f[0] = T_f_initial       # C
        results.T_C_e[0] = T_e_initial       # C
        results.P_Pa[0] = P_f_initial        # Pa
        results.h_Jkg[0] = FluidState.getHFromPT(results.P_Pa[0], results.T_C_f[0], self.fluid)
        results.rho_kgm3[0] = FluidState.getRhoFromPT(results.P_Pa[0], results.T_C_f[0], self.fluid)

        # Calculate the Friction Factor
        # Use limit of Colebrook-white equation for large Re
        ff = 0.25 * (1/np.log10(self.epsilon/(self.wellRadius * 2)/3.7))**2

        alpha_rock = self.params.k_rock/self.params.density_rock/self.params.C_rock  #D rock
        t_d = alpha_rock*time_seconds/(self.wellRadius**2)  #dim
        if t_d < 2.8:
            beta = ((np.pi*t_d)**-0.5 + 0.5 - 0.25*(t_d/np.pi)**0.5 + 0.125*t_d)
        else:
            beta = (2/(np.log(4*t_d)-2*0.58) - 2*0.58/(np.log(4*t_d)-2*0.58)**2)

        # loop over all well segments
        for i in range(1, self.params.N_dx+1):
            results.z_m[i] = results.z_m[i-1] + dz

            # far-field rock temp
            results.T_C_e[i] = results.T_C_e[0] - results.z_m[i]*self.dT_dz
            # fluid velocity
            results.v_ms[i] = m_dot / A_c / results.rho_kgm3[i-1]  #m/s

            # Calculate Pressure
            results.delta_P_loss[i] =  ff * dL / ( 2 * self.wellRadius) * \
                                    results.rho_kgm3[i-1] * results.v_ms[i]**2. / 2.  #Pa
            results.rho_kgm3[i] = FluidState.getRhoFromPh(results.P_Pa[i-1], results.h_Jkg[i-1], self.fluid)
            results.P_Pa[i]         = results.P_Pa[i-1] - results.rho_kgm3[i]* self.params.g * dz \
                                    - results.delta_P_loss[i]
            # If pressure is very near and below the critical point of 7.377 MPa (for CO2),
            # there will be a coolprop convergence error
            if self.fluid == 'CO2' and results.P_Pa[i] < 7.38e6 and results.P_Pa[i] > 7.37e6:
                print('Semi_analytic_well: Manually adjusting pressure from %s' \
                        ' MPa to 7.37 MPa to avoid CoolProp CO2 critical point' \
                        ' convergence issues.'%(results.P_Pa[i]/1e6))
                results.P_Pa[i] = 7.37e6

            # Throw exception if below saturation pressure of water at previous temperature
            if self.fluid == 'Water':
                P_sat = FluidState.getPFromTQ(results.T_C_f[i-1], 0, self.fluid)
                if results.P_Pa[i] < P_sat:
                    raise ValueError('SemiAnalyticalWell:BelowSaturationPressure - ',\
                    'Below saturation pressure of water at %s m !' %(results.z_m[i]))

            h_noHX = results.h_Jkg[i-1] - self.params.g * dz
            T_noHX = FluidState.getTFromPh(results.P_Pa[i], h_noHX, self.fluid)
            results.cp_JK[i] = FluidState.getCpFromPh(results.P_Pa[i], h_noHX, self.fluid)

            #Find Fluid Temp
            if not self.params.useWellboreHeatLoss:
                results.T_C_f[i] = T_noHX
                results.h_Jkg[i] = h_noHX
                results.q[i] = 0.
            else:
                # See Zhang, Pan, Pruess, Finsterle (2011). A time-convolution
                # approach for modeling heat exchange between a wellbore and
                # surrounding formation. Geothermics 40, 261-266.
                x = dL * P_c * self.params.k_rock * beta / self.wellRadius
                y = m_dot * results.cp_JK[i]
                if math.isinf(x):
                    results.T_C_f[i] = results.T_C_e[i]
                else:
                    results.T_C_f[i] = (y * T_noHX + x *results. T_C_e[i]) / (x + y)
                results.q[i] = y * (T_noHX - results.T_C_f[i])
                results.h_Jkg[i] = FluidState.getHFromPT(results.P_Pa[i], results.T_C_f[i], self.fluid)
        return results
