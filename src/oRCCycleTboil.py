import os
import math
import numpy as np

from utils.globalConstants import getProjectRoot
from utils.fluidStates import FluidState
from src.parasiticPowerFractionCoolingTower import parasiticPowerFractionCoolingTower
from src.oRCCycleResults import ORCCycleResults

class ORCCycleTboilOutput(object):
    """ORCCycleTboilOutput."""
    pass

class ORCCycleTboil(object):
    """ ORCCycleTboil.
    Heat/power is output as specific heat and specific work. To find the
    actual power, multiply by the flowrate of geofluid through the system."""

    def __init__(self, T_ambient_C, dT_approach, dT_pinch, eta_pump, eta_turbine, coolingMode, orcFluid):

        self.T_ambient_C = T_ambient_C
        self.dT_approach = dT_approach
        self.dT_pinch = dT_pinch
        self.eta_pump = eta_pump
        self.eta_turbine = eta_turbine
        self.coolingMode = coolingMode
        self.orcFluid = orcFluid
        self.filepath = getTboilOptimum(orcFluid)

    def solve(self, T_in_C, T_boil_C = False):

        results = ORCCycleResults()

        if not T_boil_C:
            data = np.genfromtxt(self.filepath, delimiter=',')
            T_boil_C = np.interp(T_in_C, data[:,0], data[:,1])
        if T_boil_C > FluidState.getTcrit(self.orcFluid):
            raise ValueError('ORC_Cycle_Tboil:Tboil_Too_Large - Boiling temperature above critical point')

        T_condense_C = self.T_ambient_C + self.dT_approach

        # create empty list to compute cycle of 6 states
        state   = [None] * 6

        #State 1 (Condenser -> Pump)
        #saturated liquid
        state[0] = FluidState.getStateFromTQ(T_condense_C, 0, self.orcFluid)

        #State 6 (Desuperheater -> Condenser)
        #saturated vapor
        state[5] = FluidState.getStateFromTQ(state[0].T_C, 1, self.orcFluid)
        # state[5].P_Pa = state[0].P_Pa

        #State 3 (Preheater -> Boiler)
        #saturated liquid
        state[2] = FluidState.getStateFromTQ(T_boil_C, 0, self.orcFluid)

        #State 4 (Boiler -> Turbine)
        #saturated vapor
        state[3] = FluidState.getStateFromTQ(state[2].T_C, 1, self.orcFluid)
        # state[3].P_Pa = state[2].P_Pa

        #State 5 (Turbine -> Desuperheater)
        h_5s = FluidState.getHFromPS(state[0].P_Pa, state[3].S_JK, self.orcFluid)
        h_5 = state[3].h_Jkg - self.eta_turbine * (state[3].h_Jkg - h_5s)
        state[4] = FluidState.getStateFromPh(state[0].P_Pa, h_5, self.orcFluid)

        # #State 2 (Pump -> Preheater)
        h_2s = FluidState.getHFromPS(state[2].P_Pa, state[0].S_JK, self.orcFluid)
        h_2 = state[0].h_Jkg - ((state[0].h_Jkg - h_2s) / self.eta_pump)
        state[1] = FluidState.getStateFromPh(state[2].P_Pa, h_2, self.orcFluid)

        # #Calculate orc heat/work
        w_pump_orc = state[0].h_Jkg - state[1].h_Jkg
        q_preheater_orc = -1 * (state[1].h_Jkg - state[2].h_Jkg)
        q_boiler_orc = -1 * (state[2].h_Jkg - state[3].h_Jkg)
        w_turbine_orc = state[3].h_Jkg - state[4].h_Jkg
        q_desuperheater_orc = -1 * (state[4].h_Jkg - state[5].h_Jkg)
        q_condenser_orc = -1 * (state[5].h_Jkg - state[0].h_Jkg)

        # Cooling Tower Parasitic load
        dT_range = state[4].T_C - state[5].T_C
        parasiticPowerFraction = parasiticPowerFractionCoolingTower(self.T_ambient_C, self.dT_approach, dT_range, self.coolingMode)
        w_cooler_orc = q_desuperheater_orc * parasiticPowerFraction('cooling')
        w_condenser_orc = q_condenser_orc * parasiticPowerFraction('condensing')

        #water (assume pressure 100 kPa above saturation)
        P_sat = FluidState.getPFromTQ(T_in_C, 0, 'Water')
        cp = FluidState.getCpFromPT(P_sat + 100e3, T_in_C, 'Water')
        #Water state 11, inlet, 12, mid, 13 exit
        T_C_11 = T_in_C;
        T_C_12 = T_boil_C + self.dT_pinch;
        #mdot_ratio = mdot_orc / mdot_water
        mdot_ratio = cp * (T_C_11 - T_C_12) / q_boiler_orc;
        T_C_13 = T_C_12 - mdot_ratio * q_preheater_orc / cp;

        # check that T_C(13) isn't below pinch constraint
        if T_C_13 < (state[1].T_C + self.dT_pinch):
            # pinch constraint is here, not at 12
            # outlet is pump temp plus pinch
            T_C_13 = state[1].T_C + self.dT_pinch
            R = q_boiler_orc / (q_boiler_orc + q_preheater_orc)
            T_C_12 = T_C_11 - (T_C_11 - T_C_13) * R
            mdot_ratio = cp * (T_C_11 - T_C_12) / q_boiler_orc

        #Calculate water heat/work
        results.w_pump = mdot_ratio * w_pump_orc
        results.q_preheater = mdot_ratio * q_preheater_orc
        results.q_boiler = mdot_ratio * q_boiler_orc
        results.w_turbine = mdot_ratio * w_turbine_orc
        results.q_desuperheater = mdot_ratio * q_desuperheater_orc
        results.q_condenser = mdot_ratio * q_condenser_orc
        results.w_cooler = mdot_ratio * w_cooler_orc
        results.w_condenser = mdot_ratio * w_condenser_orc

        results.w_net = results.w_turbine + results.w_pump + results.w_cooler + results.w_condenser

        # Calculate temperatures
        results.dT_range_CT = state[4].T_C - state[5].T_C
        dT_A_p = T_C_13 - state[1].T_C
        dT_B_p = T_C_12 - state[2].T_C
        if dT_A_p == dT_B_p:
            results.dT_LMTD_preheater = dT_A_p
        else:
            div = dT_A_p/dT_B_p
            results.dT_LMTD_preheater = (dT_A_p - dT_B_p) / (math.log(abs(div)) * np.sign(div))

        dT_A_b = T_C_12 - state[2].T_C
        dT_B_b = T_C_11 - state[3].T_C
        if dT_A_b == dT_B_b:
            results.dT_LMTD_boiler = dT_A_b
        else:
            div = dT_A_b / dT_B_b
            results.dT_LMTD_boiler = (dT_A_b - dT_B_b) / (math.log(abs(div)) * np.sign(div))

        results.end_T_C = T_C_13

        return results

    def gatherOutput(self):
        output = ORCCycleTboilOutput()
        return output
