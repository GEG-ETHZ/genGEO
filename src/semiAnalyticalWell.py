import math
import numpy as np

from src.semiAnalyticalWellResults import SemiAnalyticalWellResults

from utils.constantsAndPaths import ConversionConstants
from utils.fluidStateFromPT import FluidStateFromPT
from utils.frictionFactor import frictionFactor

class SemiAnalyticalWell(object):
    """SemiAnalyticalWell to compute heat transport with fluid flow in a well
        and analytical conduction in the surrounding rock."""

    def __init__(self, params, T_e_initial, dz_total = 0., dr_total = 0.):
        self.params = params
        # dimensions
        # dz_total is the change in elevation, negative if going into the earth
        self.dz_total = dz_total
        # dr_total is the deviation from vertical, always positive
        self.dr_total = dr_total
        # T_e_initial is the beginning rock temp
        self.T_e_initial = T_e_initial

    def solve(self, initial_state, m_dot):
        # results
        self.results = SemiAnalyticalWellResults(self.params.N_dx, self.params.working_fluid)

        m_dot = m_dot * self.params.well_multiplier

        P_f_initial = initial_state.P_Pa()
        T_f_initial = initial_state.T_C()
        time_seconds = self.params.time_years * ConversionConstants.secPerYear

        # set geometry
        dz = self.dz_total/self.params.N_dx             # m
        dr = self.dr_total/self.params.N_dx             # m
        dL = (dz**2 + dr**2)**0.5                       # m
        A_c = np.pi * self.params.well_radius**2        # m**2
        P_c = np.pi * 2 * self.params.well_radius       # m

        # set states
        self.results.T_C_f[0] = T_f_initial             # C
        self.results.T_C_e[0] = self.T_e_initial        # C
        self.results.P_Pa[0] = P_f_initial              # Pa
        self.results.h_Jkg[0] = FluidStateFromPT.getHFromPT(self.results.P_Pa[0], self.results.T_C_f[0], self.params.working_fluid)
        self.results.rho_kgm3[0] = FluidStateFromPT.getRhoFromPT(self.results.P_Pa[0], self.results.T_C_f[0], self.params.working_fluid)
        self.results.v_ms[0] = m_dot / A_c / self.results.rho_kgm3[0]  #m/s

        # Calculate the Friction Factor
        # Use Colebrook-white equation for wellbore friction loss.
        # Calculate for first element and assume constant in remainder
        ff = frictionFactor(self.params.well_radius, self.results.P_Pa[0], self.results.h_Jkg[0], \
                            m_dot, self.params.working_fluid, self.params.epsilon)

        alpha_rock = self.params.k_rock/self.params.rho_rock/self.params.c_rock  #D rock
        t_d = alpha_rock*time_seconds/(self.params.well_radius**2)  #dim
        if t_d < 2.8:
            beta = ((np.pi*t_d)**-0.5 + 0.5 - 0.25*(t_d/np.pi)**0.5 + 0.125*t_d)
        else:
            beta = (2/(np.log(4*t_d)-2*0.58) - 2*0.58/(np.log(4*t_d)-2*0.58)**2)

        # loop over all well segments
        for i in range(1, self.params.N_dx+1):
            self.results.z_m[i] = self.results.z_m[i-1] + dz

            # far-field rock temp
            self.results.T_C_e[i] = self.results.T_C_e[0] - self.results.z_m[i] * self.params.dT_dz
            # fluid velocity
            self.results.v_ms[i] = m_dot / A_c / self.results.rho_kgm3[i-1]  #m/s

            # Calculate Pressure
            self.results.delta_P_loss[i] = ff * dL / ( 2 * self.params.well_radius) * \
                                            self.results.rho_kgm3[i-1] * \
                                            self.results.v_ms[i]**2. / 2.  #Pa
            self.results.rho_kgm3[i] = FluidStateFromPT.getRhoFromPh(self.results.P_Pa[i-1], self.results.h_Jkg[i-1], self.params.working_fluid)
            self.results.P_Pa[i]     = self.results.P_Pa[i-1] - self.results.rho_kgm3[i] * \
                                        self.params.g * dz - self.results.delta_P_loss[i]

            # Throw exception if below saturation pressure of water at previous temperature
            if self.params.working_fluid.lower() == 'water':
                P_sat = FluidStateFromPT.getPFromTQ(self.results.T_C_f[i-1], 0, self.params.working_fluid)
                if self.results.P_Pa[i] < P_sat:
                    raise ValueError('SemiAnalyticalWell:BelowSaturationPressure - ' \
                    'Below saturation pressure of water at %s m !' %(self.results.z_m[i]))

            h_noHX = self.results.h_Jkg[i-1] - self.params.g * dz
            T_noHX = FluidStateFromPT.getTFromPh(self.results.P_Pa[i], h_noHX, self.params.working_fluid)
            self.results.cp_JK[i] = FluidStateFromPT.getCpFromPh(self.results.P_Pa[i], h_noHX, self.params.working_fluid)

            #Find Fluid Temp
            if not self.params.useWellboreHeatLoss:
                self.results.T_C_f[i] = T_noHX
                self.results.h_Jkg[i] = h_noHX
                self.results.q[i] = 0.
            else:
                # See Zhang, Pan, Pruess, Finsterle (2011). A time-convolution
                # approach for modeling heat exchange between a wellbore and
                # surrounding formation. Geothermics 40, 261-266.
                x = dL * P_c * self.params.k_rock * beta / self.params.well_radius
                y = m_dot * self.results.cp_JK[i]
                if math.isinf(x):
                    self.results.T_C_f[i] = self.results.T_C_e[i]
                else:
                    self.results.T_C_f[i] = (y * T_noHX + x *self.results. T_C_e[i]) / (x + y)
                self.results.q[i] = y * (T_noHX - self.results.T_C_f[i])
                self.results.h_Jkg[i] = FluidStateFromPT.getHFromPT(self.results.P_Pa[i], self.results.T_C_f[i], self.params.working_fluid)
        return FluidStateFromPT(self.results.P_Pa[i], self.results.T_C_f[i], self.params.working_fluid)

    def gatherOutput(self):
        return self.results
