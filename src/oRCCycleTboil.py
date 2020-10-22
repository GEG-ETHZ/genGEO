import os, math
import numpy as np

from src.parasiticPowerFractionCoolingTower import parasiticPowerFractionCoolingTower
from src.powerPlantOutput import PowerPlantOutput

from utils.constantsAndPaths import getTboilOptimum
from utils.fluidStates import FluidState
from utils.fluidStateFromTQ import FluidStateFromTQ
from utils.fluidStateFromPh import FluidStateFromPh
from utils.fluidStateFromPT import FluidStateFromPT

class ORCCycleTboil(object):
    """ ORCCycleTboil.
    Heat/power is output as specific heat and specific work. To find the
    actual power, multiply by the flowrate of geofluid through the system."""

    def __init__(self, params):

        self.params = params
        self.data = np.genfromtxt(getTboilOptimum(params.orc_fluid), delimiter=',')

    def solve(self, initialState, T_boil_C = False):

        T_in_C = initialState.T_C()

        if not T_boil_C:
            T_boil_C = np.interp(T_in_C, self.data[:,0], self.data[:,1])
        if T_boil_C > FluidState.getTcrit(self.params.orc_fluid):
            raise ValueError('ORC_Cycle_Tboil:Tboil_Too_Large - Boiling temperature above critical point')

        T_condense_C = self.params.T_ambient_C + self.params.dT_approach

        # create empty list to compute cycle of 6 states
        state   = [None] * 6

        #State 1 (Condenser -> Pump)
        #saturated liquid
        state[0] = FluidStateFromTQ(T_condense_C, 0, self.params.orc_fluid)

        #State 6 (Desuperheater -> Condenser)
        #saturated vapor
        state[5] = FluidStateFromTQ(state[0].T_C(), 1, self.params.orc_fluid)
        # state[5].P_Pa = state[0].P_Pa

        #State 3 (Preheater -> Boiler)
        #saturated liquid
        state[2] = FluidStateFromTQ(T_boil_C, 0, self.params.orc_fluid)

        #State 4 (Boiler -> Turbine)
        #saturated vapor
        state[3] = FluidStateFromTQ(state[2].T_C(), 1, self.params.orc_fluid)
        # state[3].P_Pa = state[2].P_Pa

        #State 5 (Turbine -> Desuperheater)
        h_5s = FluidState.getHFromPS(state[0].P_Pa(), state[3].S_JK(), self.params.orc_fluid)
        h_5 = state[3].h_Jkg() - self.params.eta_turbine_orc * (state[3].h_Jkg() - h_5s)
        state[4] = FluidStateFromPh(state[0].P_Pa(), h_5, self.params.orc_fluid)

        # #State 2 (Pump -> Preheater)
        h_2s = FluidState.getHFromPS(state[2].P_Pa(), state[0].S_JK(), self.params.orc_fluid)
        h_2 = state[0].h_Jkg() - ((state[0].h_Jkg() - h_2s) / self.params.eta_pump_orc)
        state[1] = FluidStateFromPh(state[2].P_Pa(), h_2, self.params.orc_fluid)

        # #Calculate orc heat/work
        self.w_pump_orc = state[0].h_Jkg() - state[1].h_Jkg()
        self.q_preheater_orc = -1 * (state[1].h_Jkg() - state[2].h_Jkg())
        self.q_boiler_orc = -1 * (state[2].h_Jkg() - state[3].h_Jkg())
        self.w_turbine_orc = state[3].h_Jkg() - state[4].h_Jkg()
        self.q_desuperheater_orc = -1 * (state[4].h_Jkg() - state[5].h_Jkg())
        self.q_condenser_orc = -1 * (state[5].h_Jkg() - state[0].h_Jkg())

        # Cooling Tower Parasitic load
        dT_range = state[4].T_C() - state[5].T_C()
        parasiticPowerFraction = parasiticPowerFractionCoolingTower(self.params.T_ambient_C, self.params.dT_approach, dT_range, self.params.cooling_mode)
        self.w_cooler_orc = self.q_desuperheater_orc * parasiticPowerFraction('cooling')
        self.w_condenser_orc = self.q_condenser_orc * parasiticPowerFraction('condensing')

        #water (assume pressure 100 kPa above saturation)
        P_sat = FluidState.getPFromTQ(T_in_C, 0, 'Water')
        cp = FluidState.getCpFromPT(P_sat + 100e3, T_in_C, 'Water')
        #Water state 11, inlet, 12, mid, 13 exit
        self.T_C_11 = T_in_C;
        self.T_C_12 = T_boil_C + self.params.dT_pinch;
        #mdot_ratio = mdot_orc / mdot_water
        self.mdot_ratio = cp * (self.T_C_11 - self.T_C_12) / self.q_boiler_orc;
        self.T_C_13 = self.T_C_12 - self.mdot_ratio * self.q_preheater_orc / cp;

        # check that T_C(13) isn't below pinch constraint
        if self.T_C_13 < (state[1].T_C() + self.params.dT_pinch):
            # pinch constraint is here, not at 12
            # outlet is pump temp plus pinch
            self.T_C_13 = state[1].T_C() + self.params.dT_pinch
            R = self.q_boiler_orc / (self.q_boiler_orc + self.q_preheater_orc)
            self.T_C_12 = self.T_C_11 - (self.T_C_11 - self.T_C_13) * R
            self.mdot_ratio = cp * (self.T_C_11 - self.T_C_12) / self.q_boiler_orc

        # pass states to self variable
        self.state = state

        # return temperature
        self.state_out = FluidStateFromPT(initialState.P_Pa(), self.T_C_13, self.params.working_fluid)

        return self.state_out

    def gatherOutput(self):

        output = PowerPlantOutput()

        output.state_out = self.state_out

        #Calculate water heat/work
        output.q_preheater = self.mdot_ratio * self.q_preheater_orc
        output.q_boiler = self.mdot_ratio * self.q_boiler_orc
        output.q_desuperheater = self.mdot_ratio * self.q_desuperheater_orc
        output.q_condenser = self.mdot_ratio * self.q_condenser_orc
        output.w_turbine = self.mdot_ratio * self.w_turbine_orc
        output.w_pump = self.mdot_ratio * self.w_pump_orc
        output.w_cooler = self.mdot_ratio * self.w_cooler_orc
        output.w_condenser = self.mdot_ratio * self.w_condenser_orc
        output.w_net = output.w_turbine + output.w_pump + output.w_cooler + output.w_condenser

        # Calculate temperatures
        output.dT_range_CT = self.state[4].T_C() - self.state[5].T_C()
        dT_A_p = self.T_C_13 - self.state[1].T_C()
        dT_B_p = self.T_C_12 - self.state[2].T_C()
        if dT_A_p == dT_B_p:
            output.dT_LMTD_preheater = dT_A_p
        else:
            div = dT_A_p/dT_B_p
            output.dT_LMTD_preheater = (dT_A_p - dT_B_p) / (math.log(abs(div)) * np.sign(div))

        dT_A_b = self.T_C_12 - self.state[2].T_C()
        dT_B_b = self.T_C_11 - self.state[3].T_C()
        if dT_A_b == dT_B_b:
            output.dT_LMTD_boiler = dT_A_b
        else:
            div = dT_A_b / dT_B_b
            output.dT_LMTD_boiler = (dT_A_b - dT_B_b) / (math.log(abs(div)) * np.sign(div))

        return output
