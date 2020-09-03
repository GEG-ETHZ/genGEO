import math
import numpy as np

from CoolProp.CoolProp import PropsSI

class SemiAnalyticalWell(object):
    """SemiAnalyticalWell to compute heat transport with fluid flow in a well
        and analytical conduction in the surrounding rock."""

    def __init__(self, params):
        self.params = params
        self.alpha_rock = self.params.k_rock/self.params.density_rock/self.params.C_rock  #D rock

    # INPUT functions
    def initializeWellDims(self, dz_total, dr_total, wellRadius, N_dx = 100):
        self.wellRadius = wellRadius
        self.N_dx = N_dx
        self.dz_total = dz_total
        self.dr_total = dr_total
        self.dz = self.dz_total/self.N_dx           # m
        self.dr = self.dr_total/self.N_dx           # m
        self.dL = (self.dz**2 + self.dr**2)**0.5    # m
        self.A_c = np.pi * self.wellRadius**2       # m**2
        self.P_c = np.pi * 2 * self.wellRadius      # m

    def initializeStates(self, fluid, P_f_initial, T_f_initial, T_e_initial):
        self.fluid = fluid
        self.T_f_initial = T_f_initial
        self.T_e_initial = T_e_initial
        self.P_f_initial = P_f_initial

    def initializeModel(self, epsilon, time_seconds, m_dot, dT_dz, useWellboreHeatLoss = True):
        self.epsilon = epsilon
        self.time_seconds = time_seconds
        self.m_dot = m_dot
        self.dT_dz = dT_dz
        self.useWellboreHeatLoss =  useWellboreHeatLoss

    # OUTPUT functions
    def getState(self):
        return (self.z, self.P / 1.e6, self.T_f, self.T_e)

    def getEndPressure(self):
        return self.P[-1]

    def getEndTemperature(self):
        return self.T_f[-1]

    def getEndEnthalpy(self):
        return self.h[-1]

    def getHeat(self):
        return -1. * np.sum(self.q)

    def createEmptyArrays(self):
        # create zero arrays of the size of N_dz+1 to iterate through all n+1 well segments.
        self.z              = np.zeros(self.N_dx+1)
        self.T_f            = np.zeros(self.N_dx+1)
        self.T_e            = np.zeros(self.N_dx+1)
        self.P              = np.zeros(self.N_dx+1)
        self.h              = np.zeros(self.N_dx+1)
        self.rho_fluid      = np.zeros(self.N_dx+1)
        self.q              = np.zeros(self.N_dx+1)
        self.Cp_fluid       = np.zeros(self.N_dx+1)
        self.v              = np.zeros(self.N_dx+1)
        self.delta_P_loss   = np.zeros(self.N_dx+1)

    def setInitialConditions(self):
        self.createEmptyArrays()
        self.T_f[0] = self.T_f_initial      # C
        self.T_e[0] = self.T_e_initial      # C
        self.P[0] = self.P_f_initial        # Pa
        self.h[0] = PropsSI('HMASS', 'P', self.P[0], 'T', self.T_f[0]+273.15, self.fluid)
        self.rho_fluid[0] = PropsSI('DMASS', 'P', self.P[0], 'T', self.T_f[0]+273.15, self.fluid)

    # COMPUTE functions
    def computeSolution(self):

        # setup initial conditions
        self.setInitialConditions()

        # Calculate the Friction Factor
        # Use limit of Colebrook-white equation for large Re
        ff = 0.25 * (1/np.log10(self.epsilon/(self.wellRadius * 2)/3.7))**2

        t_d = self.alpha_rock*self.time_seconds/(self.wellRadius**2)  #dim
        if t_d < 2.8:
            beta = ((np.pi*t_d)**-0.5 + 0.5 - 0.25*(t_d/np.pi)**0.5 + 0.125*t_d)
        else:
            beta = (2/(np.log(4*t_d)-2*0.58) - 2*0.58/(np.log(4*t_d)-2*0.58)**2)

        # loop over all well segments
        for i in range(1, self.N_dx+1):
            self.z[i] = self.z[i-1] + self.dz

            # far-field rock temp
            self.T_e[i] = self.T_e[0] - self.z[i]*self.dT_dz
            # fluid velocity
            self.v[i] = self.m_dot / self.A_c / self.rho_fluid[i-1]  #m/s

            # Calculate Pressure
            self.delta_P_loss[i] =  ff * self.dL / ( 2 * self.wellRadius) * \
                                    self.rho_fluid[i-1] * self.v[i]**2. / 2.  #Pa
            self.rho_fluid[i] = PropsSI('DMASS', 'P', self.P[i-1], 'HMASS', self.h[i-1], self.fluid)
            self.P[i]         = self.P[i-1] - self.rho_fluid[i]* self.params.g * self.dz \
                                    - self.delta_P_loss[i]
            # If pressure is very near and below the critical point of 7.377 MPa (for CO2),
            # there will be a coolprop convergence error
            if self.fluid == 'CO2' and self.P[i] < 7.38e6 and self.P[i] > 7.37e6:
                print('Semi_analytic_well: Manually adjusting pressure from %s' \
                        ' MPa to 7.37 MPa to avoid CoolProp CO2 critical point' \
                        ' convergence issues.'%(self.P[i]/1e6))
                self.P[i] = 7.37e6

            # Throw exception if below saturation pressure of water at previous temperature
            if self.fluid == 'Water':
                P_sat = PropsSI('P', 'T', self.T_f[i-1] + 273.15, 'Q', 0, self.fluid)
                if self.P[i] < P_sat:
                    raise ValueError('SemiAnalyticalWell:BelowSaturationPressure - ',\
                    'Below saturation pressure of water at %s m !' %(self.z[i]))

            h_noHX = self.h[i-1] - self.params.g * self.dz
            T_noHX = PropsSI('T', 'P', self.P[i], 'HMASS', h_noHX, self.fluid) - 273.15
            self.Cp_fluid[i] = PropsSI('CPMASS', 'P', self.P[i], 'HMASS', h_noHX, self.fluid)

            #Find Fluid Temp
            if not self.useWellboreHeatLoss:
                self.T_f[i] = T_noHX
                self.h[i] = h_noHX
                self.q[i] = 0.
            else:
                # See Zhang, Pan, Pruess, Finsterle (2011). A time-convolution
                # approach for modeling heat exchange between a wellbore and
                # surrounding formation. Geothermics 40, 261-266.
                x = self.dL * self.P_c * self.params.k_rock * beta / self.wellRadius
                y = self.m_dot * self.Cp_fluid[i]
                if math.isinf(x):
                    self.T_f[i] = self.T_e[i]
                else:
                    self.T_f[i] = (y * T_noHX + x *self. T_e[i]) / (x + y)
                self.q[i] = y * (T_noHX - self.T_f[i])
                self.h[i] = PropsSI('HMASS', 'P', self.P[i], 'T', self.T_f[i] + 273.15, self.fluid)
