import unittest
import numpy as np

from src.semiAnalyticalWell import SemiAnalyticalWell
from utils.globalProperties import GlobalSimulationProperties
from utils.globalConstants import globalConstants
from tests.testAssertion import testAssert

class semiAnalyticalWellTest(unittest.TestCase):

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
                                dT_dz = -0.035)

        wellresult = well.solve(P_f_initial = 25.e6,
                                T_f_initial = 97.,
                                T_e_initial = 102.5,
                                time_seconds = 10. * globalConstants.secPerYear,
                                m_dot = 136.)

        self.assertMessages(well.fluid,
                        (wellresult.end_P_Pa(), 1.2459e6),
                        (wellresult.end_T_C(), 96.075),
                        (wellresult.end_h_Jkg(), 4.0350e5))

        # CO2
        well.fluid = 'CO2'
        wellresult = well.solve(P_f_initial = 25.e6,
                    T_f_initial = 97.,
                    T_e_initial = 102.5,
                    time_seconds = 10. * globalConstants.secPerYear,
                    m_dot = 136.)

        self.assertMessages(well.fluid,
                        (wellresult.end_P_Pa(), 1.1763e7),
                        (wellresult.end_T_C(), 57.26),
                        (wellresult.end_h_Jkg(), 3.767e5))

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
                                        dT_dz = 0.06)

        vertical_well_results = vertical_well.solve(P_f_initial = 1.e6,
                            T_f_initial = 25.,
                            T_e_initial = 15.,
                            time_seconds = 10. * globalConstants.secPerYear,
                            m_dot = 5.)

        self.assertMessages(vertical_well.fluid,
                        (vertical_well_results.end_P_Pa(), 3.533e7),
                        (vertical_well_results.end_T_C(), 67.03),
                        (vertical_well_results.end_h_Jkg(), 3.0963e5))

        horizontal_well = SemiAnalyticalWell(params = gpp,
                                        dr_total = 3000.,
                                        wellRadius = 0.279,
                                        fluid = 'water',
                                        epsilon = 55 * 1e-6,
                                        dT_dz = 0.06)

        horizontal_well_results = horizontal_well.solve(P_f_initial = vertical_well_results.end_P_Pa(),
                            T_f_initial = vertical_well_results.end_T_C(),
                            T_e_initial = 15.+ horizontal_well.dT_dz * abs(vertical_well.dz_total),
                            time_seconds = 10. * globalConstants.secPerYear,
                            m_dot = 5.)

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
                                            dT_dz = 0.06)

        vertical_well_results = vertical_well.solve(P_f_initial = 1.e6,
                            T_f_initial = 25.,
                            T_e_initial = 15.,
                            time_seconds = 10. * globalConstants.secPerYear,
                            m_dot = 5.)

        self.assertMessages(vertical_well.fluid,
                        (vertical_well_results.end_P_Pa(), 1.7245e6),
                        (vertical_well_results.end_T_C(), 156.08),
                        (vertical_well_results.end_h_Jkg(), 6.1802e5))

        horizontal_well = SemiAnalyticalWell(params = gpp,
                                            dr_total = 3000.,
                                            wellRadius = 0.279,
                                            fluid = 'CO2',
                                            epsilon = 55 * 1e-6,
                                            dT_dz = 0.06)

        horizontal_well_results = horizontal_well.solve(P_f_initial = vertical_well_results.end_P_Pa(),
                            T_f_initial = vertical_well_results.end_T_C(),
                            T_e_initial = 15. + horizontal_well.dT_dz * abs(vertical_well.dz_total),
                            time_seconds = 10. * globalConstants.secPerYear,
                            m_dot = 5.)

        self.assertMessages(horizontal_well.fluid,
                        (horizontal_well_results.end_P_Pa(), 1.7238e6),
                        (horizontal_well_results.end_T_C(), 212.746),
                        (horizontal_well_results.end_h_Jkg(), 6.755e5))
