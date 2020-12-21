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

from src.oRCCycleTboil import ORCCycleTboil
from src.coolingCondensingTower import CoolingCondensingTower
from models.coolingCondensingTowerMode import CoolingCondensingTowerMode

from utils.fluidState import FluidState
from models.simulationParameters import SimulationParameters

from tests.testAssertion import testAssert

params = SimulationParameters(orc_fluid = 'R245fa')
cycle = ORCCycleTboil(params = params)

class ORCCycleTboilTest(unittest.TestCase):

    def testParasiticPowerFraction(self):
        parasiticPowerFraction = CoolingCondensingTower.parasiticPowerFraction(15. , 7. , 25., CoolingCondensingTowerMode.Wet)
        self.assertTrue(*testAssert(parasiticPowerFraction('cooling'), 0.016025303571428565, 'Wet - cooling'))
        self.assertTrue(*testAssert(parasiticPowerFraction('condensing'), 0.02685987257142855, 'Wet - condensing'))
        self.assertRaises(Exception, parasiticPowerFraction, 'heating')
        parasiticPowerFraction = CoolingCondensingTower.parasiticPowerFraction(15. , 7. , 25., CoolingCondensingTowerMode.Dry)
        self.assertTrue(*testAssert(parasiticPowerFraction('cooling'), 0.11328571428571428, 'Dry - cooling'))
        self.assertTrue(*testAssert(parasiticPowerFraction('condensing'), 0.08842857142857143, 'Dry - condensing'))
        self.assertRaises(Exception, parasiticPowerFraction, 'heating')
        self.assertRaises(Exception, CoolingCondensingTower.parasiticPowerFraction, 15. , 7. , 25., 'Mix')

    def testORCCycleTboil(self):

        initialState = FluidState.getStateFromPT(1.e6, 150., 'water')
        results = cycle.solve(initialState, T_boil_C = 100., dT_pinch = 5.)

        self.assertTrue(*testAssert(results.state.T_C, 68.36, 'test1_temp'))
        self.assertTrue(*testAssert(results.w_net, 3.8559e4, 'test1_w_net'))
        self.assertTrue(*testAssert(results.w_turbine, 4.7773e4, 'test1_w_turbine'))
        self.assertTrue(*testAssert(results.q_preheater, 1.5778e5, 'test1_q_preheater'))
        self.assertTrue(*testAssert(results.q_boiler, 1.9380e5, 'test1_q_boiler'))

    def testORCCycleTboilFail(self):

        initialState = FluidState.getStateFromPT(1.e6, 15., 'water')
        try:
            results = cycle.solve(initialState, T_boil_C = 100., dT_pinch = 5.)
        except Exception as ex:
            self.assertTrue(str(ex).find('GenGeo::ORCCycleTboil:Tboil_Too_Large') > -1, 'test1_fail_not_found')
