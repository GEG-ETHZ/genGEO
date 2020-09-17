import math
import numpy as np

from utils.fluidStates import FluidState
from utils.globalConstants import globalConstants
from src.semiAnalyticalWellResults import SemiAnalyticalWellResults

class SemiAnalyticalWell(object):
    """SemiAnalyticalWell to compute heat transport with fluid flow in a well
        and analytical conduction in the surrounding rock."""

    def __init__(self, params, wellRadius, fluid, epsilon, dT_dz, T_e_initial, dz_total = 0., dr_total = 0.):
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
        self.T_e_initial = T_e_initial
        # data
        self.results = SemiAnalyticalWellResults(params.N_dx, fluid)

    def solve(self, initial_state, m_dot, time_years):

        P_f_initial = initial_state.P_Pa
        T_f_initial = initial_state.T_C
        time_seconds = time_years * globalConstants.secPerYear

        # set geometry
        dz = self.dz_total/self.params.N_dx           # m
        dr = self.dr_total/self.params.N_dx           # m
        dL = (dz**2 + dr**2)**0.5    # m
        A_c = np.pi * self.wellRadius**2       # m**2
        P_c = np.pi * 2 * self.wellRadius      # m

        # set states
        self.results.T_C_f[0] = T_f_initial       # C
        self.results.T_C_e[0] = self.T_e_initial       # C
        self.results.P_Pa[0] = P_f_initial        # Pa
        self.results.h_Jkg[0] = FluidState.getHFromPT(self.results.P_Pa[0], self.results.T_C_f[0], self.fluid)
        self.results.rho_kgm3[0] = FluidState.getRhoFromPT(self.results.P_Pa[0], self.results.T_C_f[0], self.fluid)

        # Calculate the Friction Factor
        # Use limit of Colebrook-white equation for large Re
        ff = 0.25 * (1/np.log10(self.epsilon/(self.wellRadius * 2)/3.7))**2

        alpha_rock = self.params.k_rock/self.params.rho_rock/self.params.C_rock  #D rock
        t_d = alpha_rock*time_seconds/(self.wellRadius**2)  #dim
        if t_d < 2.8:
            beta = ((np.pi*t_d)**-0.5 + 0.5 - 0.25*(t_d/np.pi)**0.5 + 0.125*t_d)
        else:
            beta = (2/(np.log(4*t_d)-2*0.58) - 2*0.58/(np.log(4*t_d)-2*0.58)**2)

        # loop over all well segments
        for i in range(1, self.params.N_dx+1):
            self.results.z_m[i] = self.results.z_m[i-1] + dz

            # far-field rock temp
            self.results.T_C_e[i] = self.results.T_C_e[0] - self.results.z_m[i]*self.dT_dz
            # fluid velocity
            self.results.v_ms[i] = m_dot / A_c / self.results.rho_kgm3[i-1]  #m/s

            # Calculate Pressure
            self.results.delta_P_loss[i] =  ff * dL / ( 2 * self.wellRadius) * \
                                    self.results.rho_kgm3[i-1] * self.results.v_ms[i]**2. / 2.  #Pa
            self.results.rho_kgm3[i] = FluidState.getRhoFromPh(self.results.P_Pa[i-1], self.results.h_Jkg[i-1], self.fluid)
            self.results.P_Pa[i]         = self.results.P_Pa[i-1] - self.results.rho_kgm3[i]* self.params.g * dz \
                                    - self.results.delta_P_loss[i]
            # If pressure is very near and below the critical point of 7.377 MPa (for CO2),
            # there will be a coolprop convergence error
            if self.fluid == 'CO2' and self.results.P_Pa[i] < 7.38e6 and self.results.P_Pa[i] > 7.37e6:
                print('Semi_analytic_well: Manually adjusting pressure from %s' \
                        ' MPa to 7.37 MPa to avoid CoolProp CO2 critical point' \
                        ' convergence issues.'%(self.results.P_Pa[i]/1e6))
                self.results.P_Pa[i] = 7.37e6

            # Throw exception if below saturation pressure of water at previous temperature
            if self.fluid == 'Water':
                P_sat = FluidState.getPFromTQ(self.results.T_C_f[i-1], 0, self.fluid)
                if self.results.P_Pa[i] < P_sat:
                    raise ValueError('SemiAnalyticalWell:BelowSaturationPressure - ',\
                    'Below saturation pressure of water at %s m !' %(self.results.z_m[i]))

            h_noHX = self.results.h_Jkg[i-1] - self.params.g * dz
            T_noHX = FluidState.getTFromPh(self.results.P_Pa[i], h_noHX, self.fluid)
            self.results.cp_JK[i] = FluidState.getCpFromPh(self.results.P_Pa[i], h_noHX, self.fluid)

            #Find Fluid Temp
            if not self.params.useWellboreHeatLoss:
                self.results.T_C_f[i] = T_noHX
                self.results.h_Jkg[i] = h_noHX
                self.results.q[i] = 0.
            else:
                # See Zhang, Pan, Pruess, Finsterle (2011). A time-convolution
                # approach for modeling heat exchange between a wellbore and
                # surrounding formation. Geothermics 40, 261-266.
                x = dL * P_c * self.params.k_rock * beta / self.wellRadius
                y = m_dot * self.results.cp_JK[i]
                if math.isinf(x):
                    self.results.T_C_f[i] = self.results.T_C_e[i]
                else:
                    self.results.T_C_f[i] = (y * T_noHX + x *self.results. T_C_e[i]) / (x + y)
                self.results.q[i] = y * (T_noHX - self.results.T_C_f[i])
                self.results.h_Jkg[i] = FluidState.getHFromPT(self.results.P_Pa[i], self.results.T_C_f[i], self.fluid)
        return FluidState.getStateFromPT(self.results.P_Pa[i], self.results.T_C_f[i], self.fluid)

    def gatherOutput(self):
        return self.results
