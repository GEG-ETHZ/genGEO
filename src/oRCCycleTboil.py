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
import os, math
import numpy as np

from src.coolingCondensingTower import CoolingCondensingTower
from src.powerPlantOutput import PowerPlantOutput

from utils.constantsAndPaths import getTboilOptimum
from utils.fluidState import FluidState
from utils.maxSubcritORCBoilTemp import maxSubcritORCBoilTemp
from models.simulationParameters import SimulationParameters

class ORCCycleTboil(object):
    """ ORCCycleTboil.
    Heat/power is output as specific heat and specific work. To find the
    actual power, multiply by the flowrate of geofluid through the system.
    """

    def __init__(self, params = None, **kwargs):
        self.params = params
        if self.params == None:
            self.params = SimulationParameters(**kwargs)
        self.data = getTboilOptimum()
        self.orc_fluid =  self.params.orc_fluid
        self.T_boil_max = maxSubcritORCBoilTemp(self.orc_fluid)

    def solve(self, initialState, T_boil_C = False, dT_pinch = False):

        T_in_C = initialState.T_C

        if not T_boil_C:
            T_boil_C = np.interp(T_in_C, self.data[self.params.opt_mode][self.params.orc_fluid][:,0], self.data[self.params.opt_mode][self.params.orc_fluid][:,1])

        if not dT_pinch:
            dT_pinch = np.interp(T_in_C, self.data[self.params.opt_mode][self.params.orc_fluid][:,0], self.data[self.params.opt_mode][self.params.orc_fluid][:,2])

        # run some checks  if T_in_C and T_boil_C are valid
        if np.isnan(T_in_C):
            raise Exception('GenGeo::ORCCycleTboil:T_in_NaN - ORC input temperature is NaN!')

        if np.isnan(T_boil_C):
            raise Exception('GenGeo::ORCCycleTboil:T_boil_NaN - ORC boil temperature is NaN!')

        if T_boil_C > FluidState.getTcrit(self.params.orc_fluid):
            raise Exception('GenGeo::ORCCycleTboil:Tboil_Too_Large - Boiling temperature above critical point')

        if dT_pinch <= 0:
            raise Exception('GenGeo::ORCCycleTboil:dT_pinch_Negative - dT_pinch is negative!')

        if T_in_C < T_boil_C + dT_pinch:
            raise Exception('GenGeo::ORCCycleTboil:Tboil_Too_Large - Boiling temperature of %s is greater than input temp of %s less pinch dT of %s.'%(T_boil_C, T_in_C, dT_pinch))

        # only refresh T_boil_max if orc_fluid has changed from initial
        if self.params.orc_fluid != self.orc_fluid:
            self.T_boil_max = maxSubcritORCBoilTemp(self.params.orc_fluid)
            self.orc_fluid = self.params.orc_fluid
        if T_boil_C > self.T_boil_max:
            raise Exception('GenGeo::ORCCycleTboil:Tboil_Too_Large - Boiling temperature of %s is greater than maximum allowed of %s.'%(T_boil_C, self.T_boil_max))

        T_condense_C = self.params.T_ambient_C + self.params.dT_approach

        # create empty list to compute cycle of 6 states
        state = [None] * 6

        #State 1 (Condenser -> Pump)
        #saturated liquid
        state[0] = FluidState.getStateFromTQ(T_condense_C, 0, self.params.orc_fluid)

        #State 6 (Desuperheater -> Condenser)
        #saturated vapor
        state[5] = FluidState.getStateFromTQ(state[0].T_C, 1, self.params.orc_fluid)
        # state[5].P_Pa = state[0].P_Pa

        #State 3 (Preheater -> Boiler)
        #saturated liquid
        state[2] = FluidState.getStateFromTQ(T_boil_C, 0, self.params.orc_fluid)

        #State 4 (Boiler -> Turbine)
        #saturated vapor
        state[3] = FluidState.getStateFromTQ(state[2].T_C, 1, self.params.orc_fluid)
        # state[3].P_Pa = state[2].P_Pa

        #State 5 (Turbine -> Desuperheater)
        h_5s = FluidState.getStateFromPS(state[0].P_Pa, state[3].s_JK, self.params.orc_fluid).h_Jkg
        h_5 = state[3].h_Jkg - self.params.eta_turbine_orc * (state[3].h_Jkg - h_5s)
        state[4] = FluidState.getStateFromPh(state[0].P_Pa, h_5, self.params.orc_fluid)

        # #State 2 (Pump -> Preheater)
        h_2s = FluidState.getStateFromPS(state[2].P_Pa, state[0].s_JK, self.params.orc_fluid).h_Jkg
        h_2 = state[0].h_Jkg - ((state[0].h_Jkg - h_2s) / self.params.eta_pump_orc)
        state[1] = FluidState.getStateFromPh(state[2].P_Pa, h_2, self.params.orc_fluid)

        results = PowerPlantOutput()

        # #Calculate orc heat/work
        w_pump_orc = state[0].h_Jkg - state[1].h_Jkg
        q_preheater_orc = -1 * (state[1].h_Jkg - state[2].h_Jkg)
        q_boiler_orc = -1 * (state[2].h_Jkg - state[3].h_Jkg)
        w_turbine_orc = state[3].h_Jkg - state[4].h_Jkg
        q_desuperheater_orc = -1 * (state[4].h_Jkg - state[5].h_Jkg)
        q_condenser_orc = -1 * (state[5].h_Jkg - state[0].h_Jkg)

        results.dP_pump_orc = state[1].P_Pa - state[0].P_Pa
        results.P_boil = state[2].P_Pa

        # Cooling Tower Parasitic load
        dT_range = state[4].T_C - state[5].T_C
        parasiticPowerFraction = CoolingCondensingTower.parasiticPowerFraction(self.params.T_ambient_C, self.params.dT_approach, dT_range, self.params.cooling_mode)
        w_cooler_orc = q_desuperheater_orc * parasiticPowerFraction('cooling')
        w_condenser_orc = q_condenser_orc * parasiticPowerFraction('condensing')

        #water (assume pressure 100 kPa above saturation)
        P_sat = FluidState.getStateFromTQ(T_in_C, 0, 'Water').P_Pa
        cp = FluidState.getStateFromPT(P_sat + 100e3, T_in_C, 'Water').cp_JK
        #Water state 11, inlet, 12, mid, 13 exit
        T_C_11 = T_in_C
        T_C_12 = T_boil_C + dT_pinch
        #mdot_ratio = mdot_orc / mdot_water
        mdot_ratio = cp * (T_C_11 - T_C_12) / q_boiler_orc
        T_C_13 = T_C_12 - mdot_ratio * q_preheater_orc / cp

        # check that T_C(13) isn't below pinch constraint
        if T_C_13 < (state[1].T_C + dT_pinch):
            # pinch constraint is here, not at 12
            # outlet is pump temp plus pinch
            T_C_13 = state[1].T_C + dT_pinch
            R = q_boiler_orc / (q_boiler_orc + q_preheater_orc)
            T_C_12 = T_C_11 - (T_C_11 - T_C_13) * R
            mdot_ratio = cp * (T_C_11 - T_C_12) / q_boiler_orc

        #Calculate water heat/work
        results.q_preheater = mdot_ratio * q_preheater_orc
        results.q_boiler = mdot_ratio * q_boiler_orc
        results.q_desuperheater = mdot_ratio * q_desuperheater_orc
        results.q_condenser = mdot_ratio * q_condenser_orc
        results.w_turbine = mdot_ratio * w_turbine_orc
        results.w_pump = mdot_ratio * w_pump_orc
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

        # return temperature
        results.state = FluidState.getStateFromPT(initialState.P_Pa, T_C_13, self.params.working_fluid)

        return results
