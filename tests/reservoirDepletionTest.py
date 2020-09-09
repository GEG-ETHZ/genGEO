import unittest
import numpy as np

from src.porousReservoir import PorousReservoir
from utils.depletionCurve import depletionCurve
from utils.globalProperties import GlobalSimulationProperties
from utils.globalConstants import globalConstants
from tests.testAssertion import testAssert
from utils.fluidStates import FluidState


# define global methods to be used in the tests
gpp = GlobalSimulationProperties()
gpp.k_rock = 2.1        #W/m/C
gpp.rho_rock = 2300     #kg/m^3
gpp.c_rock = 920.       #J/kg/C

reservoir = PorousReservoir(
            params = gpp,
            well_spacing = 707.,
            thickness = 300,
            permeability = 50e-15,
            T_surface_rock = 15,
            depth = 2500,
            dT_dz = 0.035,
            wellRadius = 0.205,
            reservoirConfiguration = '5spot',
            fluid = 'CO2',
            modelPressureTransient = False,
            modelTemperatureDepletion = False)

initialState = FluidState.getStateFromPT(25.e6, 40., reservoir.fluid)

class reservoirDepletionTest(unittest.TestCase):

    def testDepletionCurve(self):
        Psi_1 = 2.
        p1 = -1.9376
        p2 = 1.743
        p3 = 0.182
        gamma = np.round(depletionCurve(Psi_1, p1, p2, p3), 4)

        self.assertTrue(*testAssert(gamma, 0.2531, 'testDepletionCurve'))


    def testNoTransient(self):

        results = reservoir.solve(
                        initialState = initialState,
                        m_dot = 100,
                        time_years = 30)

        self.assertTrue(*testAssert(results.end_P_Pa, 2.0351e7, 'test2_pressure'))
        self.assertTrue(*testAssert(results.end_T_C, 102.5, 'test2_temp'))
        self.assertTrue(*testAssert(results.end_h_Jkg, 4.30e5, 'test2_enthalpy'))
        self.assertTrue(*testAssert(results.psi, 1.232, 'test2_Psi'))

    def testTransientPNoT(self):
        reservoir.modelPressureTransient = False
        reservoir.modelTemperatureDepletion = True
        results = reservoir.solve(
                        initialState = initialState,
                        m_dot = 100,
                        time_years = 30)

        self.assertTrue(*testAssert(results.end_P_Pa, 2.0351e7, 'test3_pressure'))
        self.assertTrue(*testAssert(results.end_T_C, 70.78, 'test3_temp'))
        self.assertTrue(*testAssert(results.end_h_Jkg, 3.5084e5, 'test3_enthalpy'))
        self.assertTrue(*testAssert(results.psi, 0.9082, 'test3_Psi'))

    def testTransientPT(self):
        reservoir.modelPressureTransient = True
        reservoir.modelTemperatureDepletion = True
        results = reservoir.solve(
                        initialState = initialState,
                        m_dot = 100,
                        time_years = 30)

        self.assertTrue(*testAssert(results.end_P_Pa, 1.9325e7, 'test4_pressure'))
        self.assertTrue(*testAssert(results.end_T_C, 70.78, 'test4_temp'))
        self.assertTrue(*testAssert(results.end_h_Jkg, 3.5480e5, 'test4_enthalpy'))
        self.assertTrue(*testAssert(results.psi, 0.9082, 'test4_Psi'))

    def testTransientTNoP(self):
        reservoir.modelPressureTransient = False
        reservoir.modelTemperatureDepletion = True
        results = reservoir.solve(
                        initialState = initialState,
                        m_dot = 200,
                        time_years = 30)

        self.assertTrue(*testAssert(results.end_P_Pa, 1.5702e7, 'test5_pressure'))
        self.assertTrue(*testAssert(results.end_T_C, 53.66, 'test5_temp'))
        self.assertTrue(*testAssert(results.end_h_Jkg, 3.2196e5, 'test5_enthalpy'))
        self.assertTrue(*testAssert(results.psi, 1.3701, 'test5_Psi'))
