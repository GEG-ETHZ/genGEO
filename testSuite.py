import unittest

from tests.semiAnalyticalWellTest import semiAnalyticalWellTest
from tests.reservoirDepletionTest import reservoirDepletionTest
from tests.oRCCycleTboilTest import oRCCycleTboilTest
# from tests.oRCCycleSupercritPboilTest import oRCCycleSupercritPboilTest
# from tests.oRCCycleSupercritPboilTest import heatExchangerTest
from tests.fluidSystemWaterTest import FluidSystemWaterTest

# # TODO: this is not considered right now. Use runner = unittest.TextTestRunner() to use it
# # TODO: if runner = unittest.TextTestRunner(verbosity=2) is used no test report is written
# def testSuite():
#     suite = unittest.TestSuite()
#     suite.addTest(semiAnalyticalWellTest('testProductionWell'))
#     suite.addTest(semiAnalyticalWellTest('testInjectionWellWater'))
#     suite.addTest(semiAnalyticalWellTest('testInjectionWellCO2'))
#     return suite

if __name__ == '__main__':
    runner = unittest.main()
    runner.run(testSuite())
