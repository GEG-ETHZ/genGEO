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

from src.fullSystemORC import FullSystemORC
from src.fullSystemSolver import FullSystemSolver

from models.simulationParameters import SimulationParameters

from tests.testAssertion import testAssert



class FluidSystemWaterTest(unittest.TestCase):

    def testFluidSystemWaterSolverMdot1(self):
        params = SimulationParameters(working_fluid = 'water', orc_fluid = 'R245fa', capacity_factor = 0.9)
        params.m_dot_IP = 1

        output = FullSystemORC.getDefaultWaterSystem(params).solve()

        self.assertTrue(*testAssert(output.fluid_system_solver.pump.well.state.P_Pa, 7.84791e5, 'test_subsurface_solver1_pressure'))
        self.assertTrue(*testAssert(output.fluid_system_solver.pump.well.state.T_C, 81.2595, 'test_subsurface_solver1_temp'))
        self.assertTrue(*testAssert(output.energy_results.W_net, 5.2775e3, 'test_subsurface_solver1_w_net'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_brownfield, 9.1965e6, 'test_solver1_C_brownfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_greenfield, 2.2308e7, 'test_solver1_C_greenfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 0.0083879, 'test_solver1_LCOE_brownfield'))

    def testFluidSystemWaterSolverMdot40(self):
        params = SimulationParameters(working_fluid = 'water', orc_fluid = 'R245fa', capacity_factor = 0.9)
        params.m_dot_IP = 40

        output = FullSystemORC.getDefaultWaterSystem(params).solve()

        self.assertTrue(*testAssert(output.fluid_system_solver.pump.well.state.P_Pa, 7.9610e6, 'test_subsurface_solver2_pressure'))
        self.assertTrue(*testAssert(output.fluid_system_solver.pump.well.state.T_C, 100.3125, 'test_subsurface_solver2_temp'))
        self.assertTrue(*testAssert(output.energy_results.W_net, 1.4286e+05, 'test_subsurface_solver2_w_net'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_brownfield, 2.3287e7, 'test_solver2_C_brownfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_greenfield, 3.6399062e7, 'test_solver2_C_greenfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 7.8465e-4, 'test_solver2_LCOE_brownfield'))

    def testFluidSystemWaterSolverOptMdot(self):

        params = SimulationParameters(working_fluid = 'water', orc_fluid = 'R245fa', capacity_factor = 0.9)

        full_system = FullSystemORC.getDefaultWaterSystem(params)
        full_system_solver = FullSystemSolver(full_system)

        output = full_system_solver.solve()

        self.assertTrue(*testAssert(output.optMdot, 21.9823, 'test_optMdot_solver_optMdot', 1e-3))
        self.assertTrue(*testAssert(output.energy_results.W_net, 1.5212e5, 'test_optMdot_solver_w_net', 1e-3))
        self.assertTrue(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 5.4633e-4, 'test_optMdot_solver_brownfield', 1e-3))
