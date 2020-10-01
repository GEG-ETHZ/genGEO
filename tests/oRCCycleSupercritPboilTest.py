import unittest

from src.oRCCycleSupercritPboil import ORCCycleSupercritPboil

from utils.fluidStateFromPT import FluidStateFromPT

from tests.testAssertion import testAssert


class ORCCycleSupercritPboilTest(unittest.TestCase):
    def testORCCycleSupercritPboil(self):

        cycle = ORCCycleSupercritPboil(T_ambient_C = 15.,
                                        dT_approach = 7.,
                                        dT_pinch = 5.,
                                        eta_pump = 0.9,
                                        eta_turbine = 0.8,
                                        coolingMode = 'Wet',
                                        orcFluid = 'R245fa')

        initialState = FluidStateFromPT(1.e6, 190., 'CO2')
        results = cycle.solve(initialState = initialState,
                                P_boil_Pa = 5e6)

        output = cycle.gatherOutput()
        self.assertTrue(*testAssert(output.dT_range_CT, 6.8948, 'test1_dT_range_CT'))
        self.assertTrue(*testAssert(output.w_pump, -1.9676e+03, 'test1_w_pump'))
        self.assertTrue(*testAssert(output.q_boiler, 1.2030e+05, 'test1_q_boiler'))
        self.assertTrue(*testAssert(output.w_turbine, 2.4353e+04, 'test1_w_turbine'))
        self.assertTrue(*testAssert(output.q_recuperator, 9.8023e+03, 'test1_q_recuperator'))
        self.assertTrue(*testAssert(output.q_desuperheater, -3.0484e+03, 'test1_q_desuperheater'))
        self.assertTrue(*testAssert(output.q_condenser, -9.4878e+04, 'test1_q_condenser'))
        self.assertTrue(*testAssert(output.w_cooler, -51.2720, 'test1_w_cooler'))
        self.assertTrue(*testAssert(output.w_condenser, -2.5484e+03, 'test1_w_condenser'))
        self.assertTrue(*testAssert(output.w_net, 1.9786e+04, 'test1_w_net'))
        self.assertTrue(*testAssert(output.state_out.T_C(), 73.7974, 'test1_end_T_C'))
        self.assertTrue(*testAssert(output.dT_LMTD_boiler, 9.6698, 'test1_dT_LMTD_boiler'))
        self.assertTrue(*testAssert(output.dT_LMTD_recuperator, 7.4340, 'test1_dT_LMTD_recuperator'))
