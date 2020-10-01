import unittest

from src.oRCCycleTboil import ORCCycleTboil
from src.parasiticPowerFractionCoolingTower import parasiticPowerFractionCoolingTower

from utils.fluidStateFromPT import FluidStateFromPT

from tests.testAssertion import testAssert


class ORCCycleTboilTest(unittest.TestCase):

    def testParasiticPowerFraction(self):
        parasiticPowerFraction = parasiticPowerFractionCoolingTower(15. , 7. , 25., 'Wet')
        self.assertTrue(*testAssert(parasiticPowerFraction('cooling'), 0.016025303571428565, 'Wet - cooling'))
        self.assertTrue(*testAssert(parasiticPowerFraction('condensing'), 0.02685987257142855, 'Wet - condensing'))
        self.assertRaises(ValueError, parasiticPowerFraction, 'heating')
        parasiticPowerFraction = parasiticPowerFractionCoolingTower(15. , 7. , 25., 'Dry')
        self.assertTrue(*testAssert(parasiticPowerFraction('cooling'), 0.1714285714285714, 'Dry - cooling'))
        self.assertTrue(*testAssert(parasiticPowerFraction('condensing'), 0.08842857142857143, 'Dry - condensing'))
        self.assertRaises(ValueError, parasiticPowerFraction, 'heating')
        self.assertRaises(ValueError, parasiticPowerFractionCoolingTower, 15. , 7. , 25., 'Mix')


    def testORCCycleTboil(self):

        cycle = ORCCycleTboil(T_ambient_C = 15.,
                                dT_approach = 7.,
                                dT_pinch = 5.,
                                eta_pump = 0.9,
                                eta_turbine = 0.8,
                                coolingMode = 'Wet',
                                orcFluid = 'R245fa')

        initialState = FluidStateFromPT(1.e6, 150., 'water')
        results = cycle.solve(initialState, T_boil_C = 100.)
        output = cycle.gatherOutput()
        self.assertTrue(*testAssert(output.state_out.T_C(), 68.36, 'test1_temp'))
        self.assertTrue(*testAssert(output.w_net, 3.8559e4, 'test1_w_net'))
        self.assertTrue(*testAssert(output.w_turbine, 4.7773e4, 'test1_w_turbine'))
        self.assertTrue(*testAssert(output.q_preheater, 1.5778e5, 'test1_q_preheater'))
        self.assertTrue(*testAssert(output.q_boiler, 1.9380e5, 'test1_q_boiler'))

        initialState = FluidStateFromPT(1.e6, 15., 'water')
        results = cycle.solve(initialState, T_boil_C = 100.)
        output = cycle.gatherOutput()
        self.assertTrue(*testAssert(output.state_out.T_C(), 178.2724, 'test2_temp'))
        self.assertTrue(*testAssert(output.w_net, -7.5000e4, 'test2_w_net'))
        self.assertTrue(*testAssert(output.w_turbine, -9.2921e4, 'test2_w_turbine'))
        self.assertTrue(*testAssert(output.q_preheater, -3.0690e5, 'test2_q_preheater'))
        self.assertTrue(*testAssert(output.q_boiler, -3.7696e5, 'test2_q_boiler'))
