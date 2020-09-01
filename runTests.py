import unittest

from testsuite.semiAnalyticalWellTest import runTest


def suite():
    suite = unittest.TestSuite()
    suite.addTest(runTest('test1'))
    suite.addTest(runTest('test2'))
    return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
