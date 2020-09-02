import unittest

from tests.semiAnalyticalWellTest import semiAnalyticalWellTest


def testSuite():
    suite = unittest.TestSuite()
    suite.addTest(semiAnalyticalWellTest('testProductionWell'))
    suite.addTest(semiAnalyticalWellTest('testInjectionWellWater'))
    suite.addTest(semiAnalyticalWellTest('testInjectionWellCO2'))
    return suite

if __name__ == '__main__':
    runner = unittest.main()
    runner.run(testSuite())
