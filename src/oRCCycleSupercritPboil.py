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
from src.heatExchanger import heatExchanger
from src.heatExchangerOptMdot import heatExchangerOptMdot
from src.powerPlantOutput import PowerPlantOutput

from models.simulationParameters import SimulationParameters

from utils.constantsAndPaths import getPboilOptimum
from utils.fluidState import FluidState

class ORCCycleSupercritPboil(object):
    """ ORCCycleSupercritPboil.
    Heat/power is output as specific heat and specific work. To find the
    actual power, multiply by the flowrate of geofluid through the system.
    """
    def __init__(self, params = None, **kwargs):
        self.params = params
        if self.params == None:
            self.params = SimulationParameters(**kwargs)
        self.data = getPboilOptimum()

    def solve(self, initialState, P_boil_Pa = False):

        T_in_C = initialState.T_C

        if not P_boil_Pa:
            P_boil_Pa = np.interp(T_in_C, self.data[:,0], self.data[:,1])
        # Critical point of R245fa
        # if Pboil is below critical, throw error
        if P_boil_Pa < FluidState.getPcrit(self.params.orc_fluid):
            raise Exception('GenGeo::ORCCycleSupercritPboil:lowBoilingPressure - Boiling Pressure Below Critical Pressure')
        # The line of minimum entropy to keep the fluid vapor in turbine is
        # entropy at saturated vapor at 125C. So inlet temp must provide this
        # minimum entropy.
        s_min = FluidState.getStateFromTQ(125., 1, self.params.orc_fluid).s_JK
        T_min = FluidState.getStateFromPS(P_boil_Pa, s_min, self.params.orc_fluid).T_C
        if (T_in_C - self.params.dT_pinch) < T_min:
            raise Exception('GenGeo::ORCCycleSupercritPboil:lowInletTemp - Inlet Temp below %.1f C for Supercritical Fluid'%(T_min+self.params.dT_pinch))

        T_condense_C = self.params.T_ambient_C + self.params.dT_approach

        # create empty list to compute cycle of 7 states
        state   = [None] * 7

        #State 1 (Condenser -> Pump)
        #saturated liquid
        state[0] = FluidState.getStateFromTQ(T_condense_C, 0, self.params.orc_fluid)

        # State 7 (Desuperheater -> Condenser)
        # saturated vapor
        state[6] = FluidState.getStateFromTQ(state[0].T_C, 1, self.params.orc_fluid)

        # State 2 (Pump -> Recuperator)
        h_2s = FluidState.getStateFromPS(P_boil_Pa, state[0].s_JK, self.params.orc_fluid).h_Jkg
        h2 = state[0].h_Jkg - ((state[0].h_Jkg - h_2s) / self.params.eta_pump_orc)
        state[1] = FluidState.getStateFromPh(P_boil_Pa, h2, self.params.orc_fluid)

        # water (assume pressure 100 kPa above saturation)
        P_water = FluidState.getStateFromTQ(T_in_C, 0, 'Water').P_Pa + 100e3

        # Guess orc_in fluid is state[1].T_C
        state[2] = FluidState.getStateFromPT(state[1].P_Pa, state[1].T_C, self.params.orc_fluid)
        # Water in temp is T_in_C
        T_C_11 = T_in_C
        P_4 = state[1].P_Pa

        # initialize state 2 with pressure state 2 and dummy temperature
        state[3] = FluidState.getStateFromPT(P_4, 15., self.params.orc_fluid)

        # initialize state 3 with pressure state 1 and dummy enthalpy
        state[4] =  FluidState.getStateFromPh(state[0].P_Pa, 1e4, self.params.orc_fluid)

        # initialize state 6 with pressure state 1 and dummy temperature
        state[5] = FluidState.getStateFromPT(state[0].P_Pa, 15., self.params.orc_fluid)

        results = PowerPlantOutput()

        dT = 1
        while abs(dT) >= 1:
            # State 4 (Boiler -> Turbine)
            # Input orc/geo heat exchanger
            opt_heatExchanger_results = heatExchangerOptMdot(state[2].T_C, P_4, self.params.orc_fluid, T_C_11, P_water, 'Water', self.params.dT_pinch, T_min)
            # state[3] = FluidStateFromPT(P_4, opt_heatExchanger_results.T_1_out, self.params.orc_fluid)
            state[3].T_C = opt_heatExchanger_results.T_1_out

            #State 5 (Turbine -> Recuperator)
            h_5s = FluidState.getStateFromPS(state[0].P_Pa, state[3].s_JK, self.params.orc_fluid).h_Jkg
            h_5 = state[3].h_Jkg - self.params.eta_turbine_orc * (state[3].h_Jkg - h_5s)
            state[4].h_Jkg = h_5
            # state[4] =  FluidStateFromPh(state[0].P_Pa(), h_5, self.params.orc_fluid)

            # State 3 (Recuperator -> Boiler)
            # State 6 (Recuperator -> Desuperheater)
            # Assume m_dot for each fluid is 1, then output is specific heat
            # exchange
            heatExchanger_results = heatExchanger(state[1].T_C, state[1].P_Pa, 1, self.params.orc_fluid,
                                                  state[4].T_C, state[0].P_Pa, 1, self.params.orc_fluid, self.params.dT_pinch)

            # state[2] = FluidStateFromPT(state[2].P_Pa, state[2].T_C, self.params.orc_fluid)
            # state[5] = FluidStateFromPT(state[0].P_Pa, heatExchanger_results.T_2_out, self.params.orc_fluid)
            state[5].T_C = heatExchanger_results.T_2_out

            dT = state[2].T_C - heatExchanger_results.T_1_out
            state[2].T_C = heatExchanger_results.T_1_out

        #Calculate orc heat/work
        w_pump_orc = state[0].h_Jkg - state[1].h_Jkg
        q_boiler_orc = -1 * (state[2].h_Jkg - state[3].h_Jkg)
        w_turbine_orc = state[3].h_Jkg - state[4].h_Jkg
        q_desuperheater_orc = -1 * (state[5].h_Jkg - state[6].h_Jkg)
        q_condenser_orc = -1 * (state[6].h_Jkg - state[0].h_Jkg)

        # Cooling Tower Parasitic load
        results.dT_range_CT = state[5].T_C - state[6].T_C
        parasiticPowerFraction = CoolingCondensingTower.parasiticPowerFraction(self.params.T_ambient_C, self.params.dT_approach, results.dT_range_CT, self.params.cooling_mode)
        w_cooler_orc = q_desuperheater_orc * parasiticPowerFraction('cooling')
        w_condenser_orc = q_condenser_orc * parasiticPowerFraction('condensing')

        #Calculate water heat/work
        results.w_pump          = opt_heatExchanger_results.mdot_ratio * w_pump_orc
        results.q_boiler        = opt_heatExchanger_results.mdot_ratio * q_boiler_orc
        results.w_turbine       = opt_heatExchanger_results.mdot_ratio * w_turbine_orc
        results.q_recuperator   = opt_heatExchanger_results.mdot_ratio * heatExchanger_results.Q_exchanged
        results.q_desuperheater = opt_heatExchanger_results.mdot_ratio * q_desuperheater_orc
        results.q_condenser     = opt_heatExchanger_results.mdot_ratio * q_condenser_orc
        results.w_cooler        = opt_heatExchanger_results.mdot_ratio * w_cooler_orc
        results.w_condenser     = opt_heatExchanger_results.mdot_ratio * w_condenser_orc

        results.w_net = results.w_turbine + results.w_pump + results.w_cooler + results.w_condenser

        results.end_T_C = opt_heatExchanger_results.T_2_out
        results.dT_LMTD_boiler = opt_heatExchanger_results.dT_LMTD
        results.dT_LMTD_recuperator = heatExchanger_results.dT_LMTD

        # return temperature
        results.state = FluidState.getStateFromPT(initialState.P_Pa, opt_heatExchanger_results.T_2_out, self.params.working_fluid)

        return results
