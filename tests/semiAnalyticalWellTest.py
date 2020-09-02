import unittest
import numpy as np

from src.semiAnalyticalWell import SemiAnalyticalWell
from utils.globalProperties import *

class semiAnalyticalWellTest(unittest.TestCase):

    def assertMessages(self, fluid, max_error, pressure, temperature, enthalpy):
        self.assertTrue(np.isclose([pressure[0]], [pressure[1]], rtol=max_error),\
            '%s_Pressure %.4e Pa instead of %.4e Pa'%(fluid, *pressure))
        self.assertTrue(np.isclose([temperature[0]], [temperature[1]], rtol=max_error),\
            '%s_Temp %.4e C instead of %.4e C'%(fluid, *temperature))
        self.assertTrue(np.isclose([enthalpy[0]], [enthalpy[1]], rtol=max_error),\
            '%s_Enthalpy %.4e J instead of %.4e J'%(fluid, *enthalpy))

    def testProductionWell(self):
        ###
        #  Testing SemiAnalyticalWell for vertical production well settings
        ###
        # load global physical properties
        gpp = GlobalPhysicalProperties()
        # relative error of the computed values compared to reference values provided by badams
        accepted_num_error = 1e-4
        # Water
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
                             time_seconds = 10. * secPerYear, \
                             m_dot = 136., \
                             dT_dz = -0.035)

        self.assertMessages(well.fluid,\
                        accepted_num_error,\
                        (well.getEndPressure(), 1.2459e6),\
                        (well.getEndTemperature(), 96.075),\
                        (well.getEndEnthalpy(), 4.0350e5))

        # CO2
        well.fluid = 'CO2'
        well.computeSolution(epsilon = 55 * 1e-6, \
                             time_seconds = 10. * secPerYear, \
                             m_dot = 136., \
                             dT_dz = -0.035)

        self.assertMessages(well.fluid,\
                        accepted_num_error,\
                        (well.getEndPressure(), 1.1763e7),\
                        (well.getEndTemperature(), 57.26),\
                        (well.getEndEnthalpy(), 3.767e5))

    def testInjectionWellWater(self):
        ###
        #  Testing SemiAnalyticalWell for vertical and horizontal injection well settings
        ###
        # load global physical properties
        gpp = GlobalPhysicalProperties()
        # relative error of the computed values compared to reference values provided by badams
        accepted_num_error = 1e-4

        wellRadius = 0.279
        dT_dz = 0.06
        dr_total_h = 3000.
        dz_total_v = -3500.
        time = 10. * secPerYear
        eps = 55 * 1e-6
        m_dot = 5.

        vertical_well = SemiAnalyticalWell(gpp)
        vertical_well.initializeWellDims(N_dx = 100, \
                                            dz_total = dz_total_v, \
                                            dr_total = 0., \
                                            wellRadius = wellRadius)
        vertical_well.initializeStates(fluid = 'water', \
                                        P_f_initial = 1.e6, \
                                        T_f_initial = 25., \
                                        T_e_initial = 15.)
        vertical_well.computeSolution(epsilon = eps, \
                                        time_seconds = time, \
                                        m_dot = m_dot, \
                                        dT_dz = dT_dz)

        self.assertMessages(vertical_well.fluid,\
                        accepted_num_error,\
                        (vertical_well.getEndPressure(), 3.533e7),\
                        (vertical_well.getEndTemperature(), 67.03),\
                        (vertical_well.getEndEnthalpy(), 3.0963e5))


        horizontal_well = SemiAnalyticalWell(gpp)
        horizontal_well.initializeWellDims(N_dx = 100, \
                                            dz_total = 0., \
                                            dr_total = dr_total_h, \
                                            wellRadius = wellRadius)
        horizontal_well.initializeStates(fluid = 'water', \
                                        P_f_initial = vertical_well.getEndPressure(), \
                                        T_f_initial = vertical_well.getEndTemperature(), \
                                        T_e_initial = 15. + dT_dz * abs(dz_total_v))
        horizontal_well.computeSolution(epsilon = eps, \
                                        time_seconds = time, \
                                        m_dot = m_dot, \
                                        dT_dz = dT_dz)

        self.assertMessages(horizontal_well.fluid,\
                        accepted_num_error,\
                        (horizontal_well.getEndPressure(), 3.533e7),\
                        (horizontal_well.getEndTemperature(), 121.99),\
                        (horizontal_well.getEndEnthalpy(), 5.3712e5))

    def testInjectionWellCO2(self):
        ###
        #  Testing SemiAnalyticalWell for horizontal injection well settings
        ###
        # load global physical properties
        gpp = GlobalPhysicalProperties()
        # relative error of the computed values compared to reference values provided by badams
        accepted_num_error = 1e-4

        wellRadius = 0.279
        dT_dz = 0.06
        dr_total_h = 3000.
        dz_total_v = -3500.
        time = 10. * secPerYear
        eps = 55 * 1e-6
        m_dot = 5.

        vertical_well = SemiAnalyticalWell(gpp)
        vertical_well.initializeWellDims(N_dx = 100, \
                                            dz_total = dz_total_v, \
                                            dr_total = 0., \
                                            wellRadius = wellRadius)
        vertical_well.initializeStates(fluid = 'CO2', \
                                        P_f_initial = 1.e6, \
                                        T_f_initial = 25., \
                                        T_e_initial = 15.)
        vertical_well.computeSolution(epsilon = eps, \
                                        time_seconds = time, \
                                        m_dot = m_dot, \
                                        dT_dz = dT_dz)

        self.assertMessages(vertical_well.fluid,\
                        accepted_num_error,\
                        (vertical_well.getEndPressure(), 1.7245e6),\
                        (vertical_well.getEndTemperature(), 156.08),\
                        (vertical_well.getEndEnthalpy(), 6.1802e5))


        horizontal_well = SemiAnalyticalWell(gpp)
        horizontal_well.initializeWellDims(N_dx = 100, \
                                            dz_total = 0., \
                                            dr_total = dr_total_h, \
                                            wellRadius = wellRadius)
        horizontal_well.initializeStates(fluid = 'CO2', \
                                        P_f_initial = vertical_well.getEndPressure(), \
                                        T_f_initial = vertical_well.getEndTemperature(), \
                                        T_e_initial = 15. + dT_dz * abs(dz_total_v))
        horizontal_well.computeSolution(epsilon = eps, \
                                        time_seconds = time, \
                                        m_dot = m_dot, \
                                        dT_dz =dT_dz)

        self.assertMessages(horizontal_well.fluid,\
                        accepted_num_error,\
                        (horizontal_well.getEndPressure(), 1.7238e6),\
                        (horizontal_well.getEndTemperature(), 212.746),\
                        (horizontal_well.getEndEnthalpy(), 6.755e5))
