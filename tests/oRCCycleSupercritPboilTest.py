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
import unittest

from src.oRCCycleSupercritPboil import ORCCycleSupercritPboil

from utils.fluidState import FluidState
from models.simulationParameters import SimulationParameters

from tests.testAssertion import testAssert


class ORCCycleSupercritPboilTest(unittest.TestCase):
    def testORCCycleSupercritPboil(self):

        params = SimulationParameters(orc_fluid = 'R245fa')

        cycle = ORCCycleSupercritPboil(params = params)

        initialState = FluidState.getStateFromPT(1.e6, 190., 'water')
        results = cycle.solve(initialState = initialState,
                                P_boil_Pa = 5e6)

        self.assertTrue(*testAssert(results.dT_range_CT, 6.8948, 'test1_dT_range_CT'))
        self.assertTrue(*testAssert(results.w_pump, -1.9676e+03, 'test1_w_pump'))
        self.assertTrue(*testAssert(results.q_boiler, 1.2030e+05, 'test1_q_boiler'))
        self.assertTrue(*testAssert(results.w_turbine, 2.4353e+04, 'test1_w_turbine'))
        self.assertTrue(*testAssert(results.q_recuperator, 9.8023e+03, 'test1_q_recuperator'))
        self.assertTrue(*testAssert(results.q_desuperheater, -3.0484e+03, 'test1_q_desuperheater'))
        self.assertTrue(*testAssert(results.q_condenser, -9.4878e+04, 'test1_q_condenser'))
        self.assertTrue(*testAssert(results.w_cooler, -51.2720, 'test1_w_cooler'))
        self.assertTrue(*testAssert(results.w_condenser, -2.5484e+03, 'test1_w_condenser'))
        self.assertTrue(*testAssert(results.w_net, 1.9786e+04, 'test1_w_net'))
        self.assertTrue(*testAssert(results.state.T_C, 73.7974, 'test1_end_T_C'))
        self.assertTrue(*testAssert(results.dT_LMTD_boiler, 9.6698, 'test1_dT_LMTD_boiler'))
        self.assertTrue(*testAssert(results.dT_LMTD_recuperator, 7.4340, 'test1_dT_LMTD_recuperator'))
