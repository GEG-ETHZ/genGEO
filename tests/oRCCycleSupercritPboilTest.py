import unittest

from src.oRCCycleSupercritPboil import ORCCycleSupercritPboil
from src.heatExchanger import heatExchanger
from src.heatExchangerOptMdot import heatExchangerOptMdot
from tests.testAssertion import testAssert
from utils.fluidStates import FluidState


class heatExchangerTest(unittest.TestCase):

    def testheatExchanger(self):

        T_1_in = 50
        P_1 = 4e6
        m_dot_1 = 1.5
        fluid_1 = 'R245fa'
        T_2_in = 170
        m_dot_2 = 1
        fluid_2 = 'Water'
        P_2 = FluidState.getPFromTQ(T_2_in, 0, fluid_2) + 100e3
        dT_pinch = 5

        results = heatExchanger(T_1_in, P_1, m_dot_1, fluid_1, T_2_in, P_2, m_dot_2, fluid_2, dT_pinch)

        self.assertTrue(*testAssert(results.Q_exchanged, 303373.3383, 'testHeatExchanger_Q_exchanged'))
        self.assertTrue(*testAssert(results.dT_LMTD, 11.48176, 'testHeatExchanger_dT_LMTD'))
        self.assertTrue(*testAssert(results.T_1_out, 159.2937, 'testHeatExchanger_T_1_out'))
        self.assertTrue(*testAssert(results.T_2_out, 99.01979, 'testHeatExchanger_T_2_out'))
        self.assertTrue(*testAssert(results.T_1[3], 71.79053, 'testHeatExchanger_T_1'))
        self.assertTrue(*testAssert(results.T_2[3], 109.80985, 'testHeatExchanger_T_2'))
        self.assertTrue(*testAssert(results.Q[3], 45526.17687, 'testHeatExchanger_Q'))

    # def testheatExchangerOptMdot(self):
    #
    #     T_1_in = 50.
    #     P_1 = 4.e6
    #     fluid_1 = 'R245fa'
    #     T_2_in = 170.
    #     fluid_2 = 'Water'
    #     P_2 = FluidState.getPFromTQ(T_2_in, 0, fluid_2) + 100e3
    #     dT_pinch = 5.
    #     T_min = 165.
    #     maximizeHeatFromStream = '2'
    #
    #     results = heatExchangerOptMdot(T_1_in, P_1, fluid_1, T_2_in, P_2, fluid_2, dT_pinch, T_min, maximizeHeatFromStream)
    #     print(results.Q_exchanged)
    #     print(results.dT_LMTD)
    #     print(results.T_1_out)
    #     print(results.T_2_out)
    #     print(results.T_1[3])
    #     print(results.T_2[3])
    #     print(results.Q[3])
    #     print(results.q_exchanged_1)
    #     print(results.q_exchanged_2)
    #     print(results.mdot_ratio)
    #     print(results.m_dot_1)
    #     print(results.m_dot_2)
    #     self.assertTrue(*testAssert(results.Q_exchanged, 303373.3383, 'testHeatExchangerOptMdot_Q_exchanged'))
    #     self.assertTrue(*testAssert(results.dT_LMTD, 11.48176, 'testHeatExchangerOptMdot_dT_LMTD'))
    #     self.assertTrue(*testAssert(results.T_1_out, 159.2937, 'testHeatExchangerOptMdot_T_1_out'))
    #     self.assertTrue(*testAssert(results.T_2_out, 99.01979, 'testHeatExchangerOptMdot_T_2_out'))
    #     self.assertTrue(*testAssert(results.T_1[3], 71.79053, 'testHeatExchangerOptMdot_T_1'))
    #     self.assertTrue(*testAssert(results.T_2[3], 109.80985, 'testHeatExchangerOptMdot_T_2'))
    #     self.assertTrue(*testAssert(results.Q[3], 45526.17687, 'testHeatExchangerOptMdot_Q'))
    #     self.assertTrue(*testAssert(results.q_exchanged_1, 45526, 'testHeatExchangerOptMdot_q_exchanged_1'))
    #     self.assertTrue(*testAssert(results.q_exchanged_2, 45526, 'testHeatExchangerOptMdot_q_exchanged_2'))
    #     self.assertTrue(*testAssert(results.mdot_ratio, 45526, 'testHeatExchangerOptMdot_mdot_ratio'))
    #     self.assertTrue(*testAssert(results.m_dot_1, 45526, 'testHeatExchangerOptMdot_m_dot_1'))
    #     self.assertTrue(*testAssert(results.m_dot_2, 45526, 'testHeatExchangerOptMdot_m_dot_2'))


class oRCCycleSupercritPboilTest(unittest.TestCase):
    def testORCCycleSupercritPboil(self):

        cycle = ORCCycleSupercritPboil(T_ambient_C = 15.,
                                        dT_approach = 7.,
                                        dT_pinch = 5.,
                                        eta_pump = 0.9,
                                        eta_turbine = 0.8,
                                        coolingMode = 'Wet',
                                        orcFluid = 'R245fa',
                                        P_boil_Pa = 5e6)

        results = cycle.solve(T_in_C = 190.)


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
        self.assertTrue(*testAssert(results.end_T_C, 73.7974, 'test1_end_T_C'))
        self.assertTrue(*testAssert(results.dT_LMTD_boiler, 9.6698, 'test1_dT_LMTD_boiler'))
        self.assertTrue(*testAssert(results.dT_LMTD_recuperator, 7.4340, 'test1_dT_LMTD_recuperator'))
