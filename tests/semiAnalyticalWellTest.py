import unittest
import numpy as np

from src.semiAnalyticalWell import SemiAnalyticalWell
from utils.globalProperties import GlobalSimulationProperties
from utils.globalConstants import globalConstants
from utils.fluidStateFromPT import FluidStateFromPT
from tests.testAssertion import testAssert

class SemiAnalyticalWellTest(unittest.TestCase):

    def assertMessages(self, fluid, pressure, temperature, enthalpy):
        self.assertTrue(*testAssert(*(*pressure, '%s_Pressure'%fluid)))
        self.assertTrue(*testAssert(*(*temperature, '%s_Temp'%fluid)))
        self.assertTrue(*testAssert(*(*enthalpy, '%s_Enthalpy'%fluid)))

    def testProductionWell(self):
        ###
        #  Testing SemiAnalyticalWell for vertical production well settings
        ###
        # load global physical properties
        gpp = GlobalSimulationProperties()
        # Water
        well = SemiAnalyticalWell(params = gpp,
                                dz_total = 2500.,
                                wellRadius = 0.205,
                                fluid = 'water',
                                epsilon = 55 * 1e-6,
                                dT_dz = -0.035,
                                T_e_initial = 102.5)

        # initial state
        initial_state = FluidStateFromPT(25.e6, 97., well.fluid)

        well.solve(initial_state = initial_state,
                                m_dot = 136.,
                                time_years = 10.)

        wellresult = well.gatherOutput()
        self.assertMessages(well.fluid,
                        (wellresult.end_P_Pa(), 1.2432e6),
                        (wellresult.end_T_C(), 96.076),
                        (wellresult.end_h_Jkg(), 4.0350e5))

        # CO2
        well.fluid = 'CO2'
        # initial state
        initial_state = FluidStateFromPT(25.e6, 97., well.fluid)
        well.solve(initial_state = initial_state,
                    time_years = 10.,
                    m_dot = 136.)

        wellresult = well.gatherOutput()
        self.assertMessages(well.fluid,
                        (wellresult.end_P_Pa(), 1.1762e7),
                        (wellresult.end_T_C(), 57.257),
                        (wellresult.end_h_Jkg(), 3.7672e5))

    def testInjectionWellWater(self):
        ###
        #  Testing SemiAnalyticalWell for vertical and horizontal injection well settings
        ###
        # load global physical properties
        gpp = GlobalSimulationProperties()

        vertical_well = SemiAnalyticalWell(params = gpp,
                                        dz_total = -3500.,
                                        wellRadius = 0.279,
                                        fluid = 'water',
                                        epsilon = 55 * 1e-6,
                                        dT_dz = 0.06,
                                        T_e_initial = 15.)

        # initial state
        initial_state = FluidStateFromPT(1.e6, 25., vertical_well.fluid)
        vertical_well.solve(initial_state = initial_state,
                                                    m_dot = 5.,
                                                    time_years = 10.)

        vertical_well_results = vertical_well.gatherOutput()
        self.assertMessages(vertical_well.fluid,
                        (vertical_well_results.end_P_Pa(), 3.533e7),
                        (vertical_well_results.end_T_C(), 67.03),
                        (vertical_well_results.end_h_Jkg(), 3.0963e5))

        horizontal_well = SemiAnalyticalWell(params = gpp,
                                        dr_total = 3000.,
                                        wellRadius = 0.279,
                                        fluid = 'water',
                                        epsilon = 55 * 1e-6,
                                        dT_dz = 0.06,
                                        T_e_initial = 15.+ 0.06 * abs(vertical_well.dz_total))

        # initial state
        initial_state = FluidStateFromPT(vertical_well_results.end_P_Pa(), vertical_well_results.end_T_C(), vertical_well_results.fluid)
        horizontal_well.solve(initial_state = initial_state,
                            m_dot = 5.,
                            time_years = 10.)

        horizontal_well_results = horizontal_well.gatherOutput()
        self.assertMessages(horizontal_well.fluid,
                        (horizontal_well_results.end_P_Pa(), 3.533e7),
                        (horizontal_well_results.end_T_C(), 121.99),
                        (horizontal_well_results.end_h_Jkg(), 5.3712e5))

    def testInjectionWellCO2(self):
        ###
        #  Testing SemiAnalyticalWell for horizontal injection well settings
        ###
        # load global physical properties
        gpp = GlobalSimulationProperties()

        vertical_well = SemiAnalyticalWell(params = gpp,
                                            dz_total = -3500.,
                                            wellRadius = 0.279,
                                            fluid = 'CO2',
                                            epsilon = 55 * 1e-6,
                                            dT_dz = 0.06,
                                            T_e_initial = 15.)

        # initial state
        initial_state = FluidStateFromPT(1.e6, 25., vertical_well.fluid)
        vertical_well.solve(initial_state = initial_state,
                                                    m_dot = 5.,
                                                    time_years = 10.)

        vertical_well_results = vertical_well.gatherOutput()
        self.assertMessages(vertical_well.fluid,
                        (vertical_well_results.end_P_Pa(), 1.7245e6),
                        (vertical_well_results.end_T_C(), 156.08),
                        (vertical_well_results.end_h_Jkg(), 6.1802e5))

        horizontal_well = SemiAnalyticalWell(params = gpp,
                                            dr_total = 3000.,
                                            wellRadius = 0.279,
                                            fluid = 'CO2',
                                            epsilon = 55 * 1e-6,
                                            dT_dz = 0.06,
                                            T_e_initial = 15. + 0.06 * abs(vertical_well.dz_total))

        # initial state
        initial_state = FluidStateFromPT(vertical_well_results.end_P_Pa(), vertical_well_results.end_T_C(), vertical_well.fluid)
        horizontal_well.solve(initial_state = initial_state,
                            m_dot = 5.,
                            time_years = 10.)

        horizontal_well_results = horizontal_well.gatherOutput()
        self.assertMessages(horizontal_well.fluid,
                        (horizontal_well_results.end_P_Pa(), 1.7235e6),
                        (horizontal_well_results.end_T_C(), 212.746),
                        (horizontal_well_results.end_h_Jkg(), 6.755e5))

    def testInjectionWellCO2HighQ(self):
        ###
        #  Testing SemiAnalyticalWell for horizontal injection well settings
        ###
        # load global physical properties
        gpp = GlobalSimulationProperties()

        vertical_well = SemiAnalyticalWell(params = gpp,
                                            dz_total = -3500.,
                                            wellRadius = 0.279,
                                            fluid = 'CO2',
                                            epsilon = 55 * 1e-6,
                                            dT_dz = 0.06,
                                            T_e_initial = 15.)

        # initial state
        initial_state = FluidStateFromPT(1.e6, 25., vertical_well.fluid)
        vertical_well.solve(initial_state = initial_state,
                                                    m_dot = 150.,
                                                    time_years = 10.)

        vertical_well_results = vertical_well.gatherOutput()
        self.assertMessages(vertical_well.fluid,
                        (vertical_well_results.end_P_Pa(), 5.1170e5),
                        (vertical_well_results.end_T_C(), 63.9786),
                        (vertical_well_results.end_h_Jkg(), 5.3677e5))

    def testInjectionWellCO2SmallWellR(self):
        ###
        #  Testing SemiAnalyticalWell for horizontal injection well settings
        ###
        # load global physical properties
        gpp = GlobalSimulationProperties()

        vertical_well = SemiAnalyticalWell(params = gpp,
                                            dz_total = -3500.,
                                            wellRadius = 0.02,
                                            fluid = 'CO2',
                                            epsilon = 55 * 1e-6,
                                            dT_dz = 0.06,
                                            T_e_initial = 15.)

        # initial state
        initial_state = FluidStateFromPT(1.e6, 25., vertical_well.fluid)
        vertical_well.solve(initial_state = initial_state,
                                                    m_dot = 0.1,
                                                    time_years = 10.)

        vertical_well_results = vertical_well.gatherOutput()
        self.assertMessages(vertical_well.fluid,
                        (vertical_well_results.end_P_Pa(), 1.1431e6),
                        (vertical_well_results.end_T_C(), 222.2246),
                        (vertical_well_results.end_h_Jkg(), 6.8717e5))
