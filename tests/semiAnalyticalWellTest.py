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

from src.semiAnalyticalWell import SemiAnalyticalWell
from models.simulationParameters import SimulationParameters
from utils.fluidState import FluidState
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
        # Water
        params = SimulationParameters(working_fluid = 'water',
                                        time_years = 10.,
                                        m_dot_IP=136.)

        # initial state
        initial_state = FluidState.getStateFromPT(25.e6, 97., params.working_fluid)
        
        well = SemiAnalyticalWell(params, T_e_initial=102.5, dz_total=2500.)
        wellresult = well.solve(initial_state)

        self.assertMessages(well.params.working_fluid,
                        (wellresult.state.P_Pa, 1.2372e6),
                        (wellresult.state.T_C, 95.0133),
                        (wellresult.state.h_Jkg, 3.9902e5))

        # CO2
        well.params.working_fluid = 'CO2'
        # initial state
        wellresult = well.solve(initial_state)

        self.assertMessages(well.params.working_fluid,
                        (wellresult.state.P_Pa, 1.1674e7),
                        (wellresult.state.T_C, 55.9783),
                        (wellresult.state.h_Jkg, 3.7225e5))

    def testInjectionWellWater(self):
        ###
        #  Testing SemiAnalyticalWell for vertical and horizontal injection well settings
        ###
        # load parameters
        gpp = SimulationParameters(working_fluid = 'water')
        gpp.time_years = 10.
        gpp.dT_dz = 0.06
        gpp.well_radius = 0.279
        gpp.m_dot_IP = 5.

        vertical_well = SemiAnalyticalWell(params = gpp,
                                        dz_total = -3500.,
                                        T_e_initial = 15.)

        # initial state
        initial_state = FluidState.getStateFromPT(1.e6, 25., gpp.working_fluid)
        vertical_well_results = vertical_well.solve(initial_state)

        self.assertMessages(gpp.working_fluid,
                        (vertical_well_results.state.P_Pa, 3.533e7),
                        (vertical_well_results.state.T_C, 67.03),
                        (vertical_well_results.state.h_Jkg, 3.0963e5))

        horizontal_well = SemiAnalyticalWell(params = gpp,
                                        dr_total = 3000.,
                                        T_e_initial = vertical_well.T_e_initial + gpp.dT_dz * abs(vertical_well.dz_total))

        # initial state
        initial_state = FluidState.getStateFromPT(vertical_well_results.state.P_Pa, vertical_well_results.state.T_C, gpp.working_fluid)
        horizontal_well_results = horizontal_well.solve(initial_state)

        self.assertMessages(gpp.working_fluid,
                        (horizontal_well_results.state.P_Pa, 3.533e7),
                        (horizontal_well_results.state.T_C, 121.99),
                        (horizontal_well_results.state.h_Jkg, 5.3712e5))

    def testInjectionWellCO2(self):
        ###
        #  Testing SemiAnalyticalWell for horizontal injection well settings
        ###
        # load global physical properties
        gpp = SimulationParameters(working_fluid = 'co2')
        gpp.time_years = 10.
        gpp.dT_dz = 0.06
        gpp.well_radius = 0.279
        gpp.m_dot_IP = 5.

        vertical_well = SemiAnalyticalWell(params = gpp,
                                            dz_total = -3500.,
                                            T_e_initial = 15.)

        # initial state
        initial_state = FluidState.getStateFromPT(1.e6, 25., gpp.working_fluid)
        vertical_well_results = vertical_well.solve(initial_state)

        self.assertMessages(gpp.working_fluid,
                        (vertical_well_results.state.P_Pa, 1.7245e6),
                        (vertical_well_results.state.T_C, 156.08),
                        (vertical_well_results.state.h_Jkg, 6.1802e5))

        horizontal_well = SemiAnalyticalWell(params = gpp,
                                            dr_total = 3000.,
                                            T_e_initial = vertical_well.T_e_initial + gpp.dT_dz * abs(vertical_well.dz_total))

        # initial state
        initial_state = FluidState.getStateFromPT(vertical_well_results.state.P_Pa, vertical_well_results.state.T_C, gpp.working_fluid)
        horizontal_well_results = horizontal_well.solve(initial_state)

        self.assertMessages(gpp.working_fluid,
                        (horizontal_well_results.state.P_Pa, 1.7235e6),
                        (horizontal_well_results.state.T_C, 212.746),
                        (horizontal_well_results.state.h_Jkg, 6.755e5))

    def testInjectionWellCO2HighQ(self):
        ###
        #  Testing SemiAnalyticalWell for horizontal injection well settings
        ###
        # load global physical properties
        gpp = SimulationParameters(working_fluid = 'co2')
        gpp.time_years = 10.
        gpp.dT_dz = 0.06
        gpp.well_radius = 0.279
        gpp.m_dot_IP = 150.

        vertical_well = SemiAnalyticalWell(params = gpp,
                                            dz_total = -3500.,
                                            T_e_initial = 15.)

        # initial state
        initial_state = FluidState.getStateFromPT(1.e6, 25., gpp.working_fluid)
        vertical_well_results = vertical_well.solve(initial_state)

        self.assertMessages(gpp.working_fluid,
                        (vertical_well_results.state.P_Pa, 5.1170e5),
                        (vertical_well_results.state.T_C, 63.9786),
                        (vertical_well_results.state.h_Jkg, 5.3677e5))

    def testInjectionWellCO2SmallWellR(self):
        ###
        #  Testing SemiAnalyticalWell for horizontal injection well settings
        ###
        # load global physical properties
        gpp = SimulationParameters(working_fluid = 'co2')
        gpp.time_years = 10.
        gpp.dT_dz = 0.06
        gpp.well_radius = 0.02
        gpp.m_dot_IP = 0.1

        vertical_well = SemiAnalyticalWell(params = gpp,
                                            dz_total = -3500.,
                                            T_e_initial = 15.)

        # initial state
        initial_state = FluidState.getStateFromPT(1.e6, 25., gpp.working_fluid)
        vertical_well_results = vertical_well.solve(initial_state)

        self.assertMessages(gpp.working_fluid,
                        (vertical_well_results.state.P_Pa, 1.1431e6),
                        (vertical_well_results.state.T_C, 222.2246),
                        (vertical_well_results.state.h_Jkg, 6.8717e5))

