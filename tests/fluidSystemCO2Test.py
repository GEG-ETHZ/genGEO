import unittest
import numpy as np

from src.fullSystemCPG import FullSystemCPG
from src.fullSystemSolver import FullSystemSolver

from models.simulationParameters import SimulationParameters

from tests.testAssertion import testAssert



class FluidSystemCO2Test(unittest.TestCase):

    def testFluidSystemCO2Mdot10(self):
        params = SimulationParameters(working_fluid = 'co2', capacity_factor = 0.9)
        params.m_dot_IP = 10
        output = FullSystemCPG.getDefaultCPGSystem(params).solve()

        print(*testAssert(output.fluid_system_solver.pp.dP_surface, 5.9421e6, 'test_dP_surface'))
        print(*testAssert(output.fluid_system_solver.production_well.state.T_C, 59.2802, 'test_T_prod_surface_C'))
        print(*testAssert(output.energy_results.W_net, 9.7662e4, 'test_W_net'))
        print(*testAssert(output.capital_cost_model.C_brownfield, 1.3271e7, 'test_C_brownfield_N'))
        print(*testAssert(output.capital_cost_model.C_greenfield, 3.5901e7, 'test_C_greenfield_N'))
        print(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 6.5410e-4, 'test_LCOE_brownfield'))


    def testFluidSystemCO2Mdot80(self):
        params = SimulationParameters(working_fluid = 'co2', capacity_factor = 0.9)
        params.m_dot_IP = 80
        output = FullSystemCPG.getDefaultCPGSystem(params).solve()

        print(*testAssert(output.fluid_system_solver.pp.dP_surface, 5.6501e6, 'test_dP_surface'))
        print(*testAssert(output.fluid_system_solver.production_well.state.T_C, 59.0540, 'test_T_prod_surface_C'))
        print(*testAssert(output.energy_results.W_net, 4.9940e+05, 'test_W_net'))
        print(*testAssert(output.capital_cost_model.C_brownfield, 2.7650e7, 'test_C_brownfield_N'))
        print(*testAssert(output.capital_cost_model.C_greenfield, 5.0280e7, 'test_C_greenfield_N'))
        print(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 2.6650e-4, 'test_LCOE_brownfield'))

    def testFluidSystemCO2Mdot200(self):
        params = SimulationParameters(working_fluid = 'co2', capacity_factor = 0.9)
        params.m_dot_IP = 200
        output = FullSystemCPG.getDefaultCPGSystem(params).solve()

        print(*testAssert(output.fluid_system_solver.pp.dP_surface, 3.468765e+06, 'test_dP_surface'))
        print(*testAssert(output.fluid_system_solver.production_well.state.T_C, 47.6070, 'test_T_prod_surface_C'))
        print(*testAssert(output.energy_results.W_net, -1.3875e6, 'test_W_net'))
        print(*testAssert(output.capital_cost_model.C_brownfield, 5.0443e7, 'test_C_brownfield_N'))
        print(*testAssert(output.capital_cost_model.C_greenfield, 7.3073e7, 'test_C_greenfield_N'))
        print(np.isnan(output.capital_cost_model.LCOE_brownfield.LCOE), 'test_LCOE_brownfield')

    def testFluidSystemCO2Mdot100(self):
        params = SimulationParameters(working_fluid = 'co2', capacity_factor = 0.9)
        params.m_dot_IP = 100
        params.depth = 2400.
        params.permeability = 1e-8 / 100.
        output = FullSystemCPG.getDefaultCPGSystem(params).solve()

        print(*testAssert(output.fluid_system_solver.pp.dP_surface, 5.0019e6, 'test_dP_surface'))
        print(*testAssert(output.fluid_system_solver.production_well.state.T_C, 55.3144, 'test_T_prod_surface_C'))
        print(*testAssert(output.fluid_system_solver.pp.dP_pump, -1.0756e6, 'test_dP_pump'))
        print(*testAssert(output.energy_results.W_net, 8.2656e5, 'test_W_net'))
        print(*testAssert(output.capital_cost_model.C_brownfield, 2.6182e7, 'test_C_brownfield_N'))
        print(*testAssert(output.capital_cost_model.C_greenfield, 4.8021e7, 'test_C_greenfield_N'))
        print(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 1.5247e-4, 'test_LCOE_brownfield'))

    def testFluidSystemCO2SolverOptMdot(self):
        params = SimulationParameters(working_fluid = 'co2', capacity_factor = 0.9)
        
        full_system = FullSystemCPG.getDefaultCPGSystem(params)
        full_system_solver = FullSystemSolver(full_system)

        output = full_system_solver.solve()

        print(*testAssert(output.optMdot, 54.8493, 'test_optMdot_solver_optMdot', 1e-3))
        print(*testAssert(output.energy_results.W_net, 4.5411e5, 'test_optMdot_solver_w_net', 1e-3))
        print(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 2.3915e-4, 'test_optMdot_solver_LCOE_brownfield', 1e-3))
