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

from src.heatExchanger import heatExchanger
from src.heatExchangerOptMdot import heatExchangerOptMdot

from utils.fluidState import FluidState

from tests.testAssertion import testAssert


class HeatExchangerTest(unittest.TestCase):

    def testHeatExchanger(self):

        T_1_in = 50.
        P_1 = 4.e6
        m_dot_1 = 1.5
        fluid_1 = 'R245fa'
        T_2_in = 170.
        m_dot_2 = 1.
        fluid_2 = 'Water'
        P_2 = FluidState.getStateFromTQ(T_2_in, 0, fluid_2).P_Pa + 100e3
        dT_pinch = 5.

        results = heatExchanger(T_1_in, P_1, m_dot_1, fluid_1, T_2_in, P_2, m_dot_2, fluid_2, dT_pinch)

        self.assertTrue(*testAssert(results.Q_exchanged, 303373.3383, 'testHeatExchanger_Q_exchanged'))
        self.assertTrue(*testAssert(results.dT_LMTD, 11.48176, 'testHeatExchanger_dT_LMTD'))
        self.assertTrue(*testAssert(results.T_1_out, 159.2937, 'testHeatExchanger_T_1_out'))
        self.assertTrue(*testAssert(results.T_2_out, 99.01979, 'testHeatExchanger_T_2_out'))
        self.assertTrue(*testAssert(results.T_1[3], 71.79053, 'testHeatExchanger_T_1'))
        self.assertTrue(*testAssert(results.T_2[3], 109.80985, 'testHeatExchanger_T_2'))
        self.assertTrue(*testAssert(results.Q[3], 45526.17687, 'testHeatExchanger_Q'))

    def testHeatExchangerOptMdot(self):

        T_1_in = 50.
        P_1 = 4.e6
        fluid_1 = 'R245fa'
        T_2_in = 170.
        fluid_2 = 'Water'
        P_2 = FluidState.getStateFromTQ(T_2_in, 0, fluid_2).P_Pa + 100e3
        dT_pinch = 5.
        T_min = 165.
        maximizeHeatFromStream = '2'

        results = heatExchangerOptMdot(T_1_in, P_1, fluid_1, T_2_in, P_2, fluid_2, dT_pinch, T_min, maximizeHeatFromStream)

        self.assertTrue(*testAssert(results.Q_exchanged, 229529.519513, 'testHeatExchangerOptMdot_Q_exchanged'))
        self.assertTrue(*testAssert(results.dT_LMTD, 12.601838, 'testHeatExchangerOptMdot_dT_LMTD'))
        self.assertTrue(*testAssert(results.T_1_out, 163.475171347679, 'testHeatExchangerOptMdot_T_1_out'))
        self.assertTrue(*testAssert(results.T_2_out, 131.20944977165277, 'testHeatExchangerOptMdot_T_2_out'))
        self.assertTrue(*testAssert(results.T_1[3], 74.66059924618094, 'testHeatExchangerOptMdot_T_1'))
        self.assertTrue(*testAssert(results.T_2[3], 137.08472742745016, 'testHeatExchangerOptMdot_T_2'))
        self.assertTrue(*testAssert(results.Q[3], 34450.560624389356, 'testHeatExchangerOptMdot_Q'))
        self.assertTrue(*testAssert(results.q_exchanged_1, 229529.51951325324, 'testHeatExchangerOptMdot_q_exchanged_1'))
        self.assertTrue(*testAssert(results.q_exchanged_2, 167077.85548339653, 'testHeatExchangerOptMdot_q_exchanged_2'))
        self.assertTrue(*testAssert(results.mdot_ratio, 1.3737878, 'testHeatExchangerOptMdot_mdot_ratio'))
        self.assertTrue(*testAssert(results.m_dot_1, 1, 'testHeatExchangerOptMdot_m_dot_1'))
        self.assertTrue(*testAssert(results.m_dot_2, 1.3737878, 'testHeatExchangerOptMdot_m_dot_2'))
