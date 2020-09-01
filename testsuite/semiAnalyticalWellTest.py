import unittest

from src.semiAnalyticalWell import SemiAnalyticalWell
from utils.globalProperties import *

class runTest(unittest.TestCase):

    def check_num_error(self, eps, val, ref):
        return ((val - ref) / ref) < eps

    def assert_Messages(self, fluid, max_error, pressure, temperature, enthalpy):
        self.assertTrue(self.check_num_error(max_error, *pressure),\
            'Production_%s_Pressure %.4e Pa instead of %.4e Pa'%(fluid, *pressure))
        self.assertTrue(self.check_num_error(max_error, *temperature),\
            'Production_%s_Temp %.4e C instead of %.4e C'%(fluid, *temperature))
        self.assertTrue(self.check_num_error(max_error, *enthalpy),\
            'Production_%s_Enthalpy %.4e J instead of %.4e J'%(fluid, *enthalpy))

    def test1(self):
        # relative error of the computed values compared to reference values provided by badams
        accepted_num_error = 1e-4

        ###
        #  Testing the function for vertical production well settings
        ###
        # Water
        gpp = GlobalPhysicalProperties()
        well = SemiAnalyticalWell(gpp)
        well.initializeWellDims(N_dx = 100, \
                                dz_total = 2500., \
                                dr_total = 0., \
                                wellRadius = 0.205)
        well.initializeStates(fluid = 'water', \
                              P_f_initial = 25.e6, \
                              T_f_initial = 97., \
                              T_e_initial = 102.5)
        well.computeSolution(epsilon = 55 * 1e-6, \
                             time_seconds = 10. * (3600. * 24. * 365.), \
                             m_dot = 136., \
                             dT_dz = -0.035)

        self.assert_Messages(well.fluid,\
                        accepted_num_error,\
                        (well.getEndPressure(), 1.2459e6),\
                        (well.getEndTemperature(), 96.075),\
                        (well.getEndEnthalpy(), 4.0350e5))

        # CO2
        well.initializeStates(fluid = 'CO2', \
                              P_f_initial = 25.e6, \
                              T_f_initial = 97., \
                              T_e_initial = 102.5)
        well.computeSolution(epsilon = 55 * 1e-6, \
                             time_seconds = 10. * (3600. * 24. * 365.), \
                             m_dot = 136., \
                             dT_dz = -0.035)

        self.assert_Messages(well.fluid,\
                        accepted_num_error,\
                        (well.getEndPressure(), 1.1763e7),\
                        (well.getEndTemperature(), 57.26),\
                        (well.getEndEnthalpy(), 3.767e5))

    def test2(self):
        # relative error of the computed values compared to reference values provided by badams
        accepted_num_error = 1e-4
        ###
        #  Testing the function for vertical and horizontal injection well settings
        ###
        # Vertical well segment and water
        gpp = GlobalPhysicalProperties()
        well = SemiAnalyticalWell(gpp)
        well.initializeWellDims(N_dx = 100, \
                                dz_total = -3500., \
                                dr_total = 0., \
                                wellRadius = 0.279)
        well.initializeStates(fluid = 'water', \
                              P_f_initial = 1.e6, \
                              T_f_initial = 25., \
                              T_e_initial = 15.)
        well.computeSolution(epsilon = 55 * 1e-6, \
                             time_seconds = 10. * (3600. * 24. * 365.), \
                             m_dot = 5., \
                             dT_dz = 0.06)

        self.assert_Messages(well.fluid,\
                        accepted_num_error,\
                        (well.getEndPressure(), 3.533e7),\
                        (well.getEndTemperature(), 67.03),\
                        (well.getEndEnthalpy(), 3.0963e5))

        # Vertical well segment and CO2
        well.initializeStates(fluid = 'CO2', \
                              P_f_initial = 1.e6, \
                              T_f_initial = 25., \
                              T_e_initial = 15.)
        well.computeSolution(epsilon = 55 * 1e-6, \
                             time_seconds = 10. * (3600. * 24. * 365.), \
                             m_dot = 5., \
                             dT_dz = 0.06)

        self.assert_Messages(well.fluid,\
                        accepted_num_error,\
                        (well.getEndPressure(), 1.7245e6),\
                        (well.getEndTemperature(), 156.08),\
                        (well.getEndEnthalpy(), 6.1802e5))

        # Horizontal well segment and water
        dT_dz = 0.06                        # to compute T_e_initial for horizontal well
        dr_total = 3000.                    # to compute T_e_initial for horizontal well
        well = SemiAnalyticalWell(gpp)
        well.initializeWellDims(N_dx = 100, \
                                dz_total = 0., \
                                dr_total = dr_total, \
                                wellRadius = 0.279)
        well.initializeStates(fluid = 'water', \
                              P_f_initial = 1.e6, \
                              T_f_initial = 25., \
                              T_e_initial = 15. + dT_dz * abs(dr_total))
        well.computeSolution(epsilon = 55 * 1e-6, \
                             time_seconds = 10. * (3600. * 24. * 365.), \
                             m_dot = 5., \
                             dT_dz = dT_dz)

        self.assert_Messages(well.fluid,\
                        accepted_num_error,\
                        (well.getEndPressure(), 3.533e7),\
                        (well.getEndTemperature(), 121.99),\
                        (well.getEndEnthalpy(), 5.3712e5))

        # Horizontal well segment and CO2
        well.initializeStates(fluid = 'CO2', \
                              P_f_initial = 1.e6, \
                              T_f_initial = 25., \
                              T_e_initial = 15. + dT_dz * abs(dr_total))
        well.computeSolution(epsilon = 55 * 1e-6, \
                             time_seconds = 10. * (3600. * 24. * 365.), \
                             m_dot = 5., \
                             dT_dz = dT_dz)

        self.assert_Messages(well.fluid,\
                        accepted_num_error,\
                        (well.getEndPressure(), 1.7238e6),\
                        (well.getEndTemperature(), 212.746),\
                        (well.getEndEnthalpy(), 6.755e5))
