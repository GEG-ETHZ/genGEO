import unittest

from src.fullSystemORC import FullSystemORC

from src.fullSystemSolver import FullSystemSolverMinLCOEBrownfield

from utils.simulationParameters import SimulationParameters

from tests.testAssertion import testAssert



class FluidSystemWaterTest(unittest.TestCase):

    def testFluidSystemWaterSolverMdot1(self):
        params = SimulationParameters(working_fluid = 'water', well_multiplier = 4., orc_fluid = 'R245fa', capacity_factor = 0.9)
        params.m_dot_IP = 1
        
        output = FullSystemORC.getDefaultWaterSystem(params).solve()

        print(*testAssert(output.fluid_system_solver.pump.well.state.P_Pa(), 7.84791e5, 'test_subsurface_solver1_pressure'))
        print(*testAssert(output.fluid_system_solver.pump.well.state.T_C(), 81.2595, 'test_subsurface_solver1_temp'))
        print(*testAssert(output.energy_results.W_net, 5.2775e3, 'test_subsurface_solver1_w_net'))
        print(*testAssert(output.capital_cost_model.C_brownfield, 9.1965e6, 'test_solver1_C_brownfield_N'))
        print(*testAssert(output.capital_cost_model.C_greenfield, 2.2308e7, 'test_solver1_C_greenfield_N'))
        print(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 0.0083879, 'test_solver1_LCOE_brownfield'))

    def testFluidSystemWaterSolverMdot40(self):
        params = SimulationParameters(working_fluid = 'water', well_multiplier = 4., orc_fluid = 'R245fa', capacity_factor = 0.9)
        params.m_dot_IP = 40
        
        output = FullSystemORC.getDefaultWaterSystem(params).solve()

        print(*testAssert(output.fluid_system_solver.pump.well.state.P_Pa(), 7.9610e6, 'test_subsurface_solver2_pressure'))
        print(*testAssert(output.fluid_system_solver.pump.well.state.T_C(), 100.3125, 'test_subsurface_solver2_temp'))
        print(*testAssert(output.energy_results.W_net, 1.4286e+05, 'test_subsurface_solver2_w_net'))
        print(*testAssert(output.capital_cost_model.C_brownfield, 2.3287e7, 'test_solver2_C_brownfield_N'))
        print(*testAssert(output.capital_cost_model.C_greenfield, 3.6399062e7, 'test_solver2_C_greenfield_N'))
        print(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 7.8465e-4, 'test_solver2_LCOE_brownfield'))

    def testFluidSystemWaterSolverOptMdot(self):
        
        params = SimulationParameters(working_fluid = 'water', well_multiplier = 4., orc_fluid = 'R245fa', capacity_factor = 0.9)
        
        full_system = FullSystemORC.getDefaultWaterSystem(params)
        full_system_solver = FullSystemSolverMinLCOEBrownfield(full_system)

        output = full_system_solver.solve()

        print(*testAssert(output.optMdot, 21.9823, 'test_optMdot_solver_optMdot', 1e-3))
        print(*testAssert(output.energy_results.W_net, 1.5212e5, 'test_optMdot_solver_w_net', 1e-3))
        print(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 5.4633e-4, 'test_optMdot_solver_brownfield', 1e-3))
