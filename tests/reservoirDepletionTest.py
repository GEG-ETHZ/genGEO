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
import numpy as np

from src.porousReservoir import PorousReservoir
from utils.depletionCurve import depletionCurve
from models.simulationParameters import SimulationParameters
from tests.testAssertion import testAssert
from utils.fluidState import FluidState


reservoir = PorousReservoir(working_fluid = 'co2',
                            reservoir_thickness = 300,
                            permeability = 50e-15,
                            time_years = 30,
                            m_dot_IP = 100)

initialState = FluidState.getStateFromPT(25.e6, 40., reservoir.params.working_fluid)

class ReservoirDepletionTest(unittest.TestCase):

    def testDepletionCurve(self):
        Psi_1 = 2.
        p1 = -1.9376
        p2 = 1.743
        p3 = 0.182
        gamma = np.round(depletionCurve(Psi_1, p1, p2, p3), 4)

        self.assertTrue(*testAssert(gamma, 0.2531, 'testDepletionCurve'))


    def testNoTransient(self):
        reservoir.modelPressureTransient = False
        reservoir.modelTemperatureDepletion = False
        results = reservoir.solve(initialState = initialState)

        self.assertTrue(*testAssert(results.state.P_Pa, 2.0351e7, 'test2_pressure'))
        self.assertTrue(*testAssert(results.state.T_C, 102.5, 'test2_temp'))
        self.assertTrue(*testAssert(results.state.h_Jkg, 4.30e5, 'test2_enthalpy'))
        self.assertTrue(*testAssert(results.psi, 1.0151, 'test2_Psi'))

    def testTransientPNoT(self):
        reservoir.modelPressureTransient = False
        reservoir.modelTemperatureDepletion = True
        results = reservoir.solve(initialState = initialState)

        self.assertTrue(*testAssert(results.state.P_Pa, 2.0351e7, 'test3_pressure'))
        self.assertTrue(*testAssert(results.state.T_C, 71.7201, 'test3_temp'))
        self.assertTrue(*testAssert(results.state.h_Jkg, 3.5324e5, 'test3_enthalpy'))
        self.assertTrue(*testAssert(results.psi, 0.8457, 'test3_Psi'))

    def testTransientPT(self):
        reservoir.modelPressureTransient = True
        reservoir.modelTemperatureDepletion = True
        results = reservoir.solve(initialState = initialState)

        self.assertTrue(*testAssert(results.state.P_Pa, 1.9325e7, 'test4_pressure'))
        self.assertTrue(*testAssert(results.state.T_C, 71.7201, 'test4_temp'))
        self.assertTrue(*testAssert(results.state.h_Jkg, 3.5732e5, 'test4_enthalpy'))
        self.assertTrue(*testAssert(results.psi, 0.8457, 'test4_Psi'))

    def testTransientTNoP(self):
        reservoir.modelPressureTransient = False
        reservoir.modelTemperatureDepletion = True
        reservoir.params.m_dot_IP = 200
        results = reservoir.solve(initialState = initialState)

        self.assertTrue(*testAssert(results.state.P_Pa, 1.5702e7, 'test5_pressure'))
        self.assertTrue(*testAssert(results.state.T_C, 55.8873, 'test5_temp'))
        self.assertTrue(*testAssert(results.state.h_Jkg, 3.2876e5, 'test5_enthalpy'))
        self.assertTrue(*testAssert(results.psi, 1.2962, 'test5_Psi'))
