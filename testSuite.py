import unittest

from tests.semiAnalyticalWellTest import semiAnalyticalWellTest


# # TODO: this is not considered right now. Use runner = unittest.TextTestRunner() to use it
# # TODO: if runner = unittest.TextTestRunner() is used no test report is written
# def testSuite():
#     suite = unittest.TestSuite()
#     suite.addTest(semiAnalyticalWellTest('testProductionWell'))
#     suite.addTest(semiAnalyticalWellTest('testInjectionWellWater'))
#     suite.addTest(semiAnalyticalWellTest('testInjectionWellCO2'))
#     return suite

if __name__ == '__main__':
    runner = unittest.main()
    runner.run(testSuite())
