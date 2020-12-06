import unittest
import numpy as np

from src.semiAnalyticalWell import SemiAnalyticalWell
from src.porousReservoir import PorousReservoir
from src.subsurfaceComponents import DownHolePump
from src.oRCCycleTboil import ORCCycleTboil

from src.fluidSystemCO2 import FluidSystemCO2
from src.capitalCostSystem import CapitalCostSystem

from src.fullSystemCPG import FullSystemCPG
from src.lCOESimple import LCOESimple
from src.capitalCostSurfacePipes import CapitalCostSurfacePipes
from src.capitalCostWell import CapitalCostWell
from src.capitalCostWellField import CapitalCostWellField
from src.capitalCostExploration import CapitalCostExploration
from src.capitalCostWellStimulation import CapitalCostWellStimulation
from src.capitalCostSurfacePlantCPG import CapitalCostSurfacePlantCPG
from src.fullSystemSolver import FullSystemSolverMinLCOEBrownfield

from utils.simulationParameters import SimulationParameters
from utils.fluidStateFromPT import FluidStateFromPT

from tests.testAssertion import testAssert


# load simulaiton parameters
params = SimulationParameters(working_fluid = 'co2')
params.well_multiplier = 4.
params.capacity_factor = 0.9

fluid_system = FluidSystemCO2(params = params)
fluid_system.injection_well = SemiAnalyticalWell(params = params,
                                    dz_total = -params.depth,
                                    T_e_initial = params.T_ambient_C)
fluid_system.reservoir = PorousReservoir(params = params)
fluid_system.production_well = SemiAnalyticalWell(params = params,
                                    dz_total = params.depth,
                                    T_e_initial = params.T_ambient_C + params.dT_dz * params.depth)

capital_cost_system = CapitalCostSystem()
capital_cost_system.CapitalCost_SurfacePlant = CapitalCostSurfacePlantCPG(params = params)
capital_cost_system.CapitalCost_SurfacePipe = CapitalCostSurfacePipes.cost(params = params)
capital_cost_system.CapitalCost_Production_Well = CapitalCostWell.cO2Baseline(params = params)
capital_cost_system.CapitalCost_Injection_Well = CapitalCostWell.cO2Baseline(params = params)
capital_cost_system.CapitalCost_Wellfield = CapitalCostWellField.cO2MonitoringBaseline(params = params)
capital_cost_system.CapitalCost_Exploration = CapitalCostExploration.cO2Baseline(params = params)
capital_cost_system.CapitalCost_Stimulation = CapitalCostWellStimulation.cost()
capital_cost_system.lcoe_model = LCOESimple(params = params)

full_system = FullSystemCPG(params, fluid_system, capital_cost_system)

class FluidSystemCO2Test(unittest.TestCase):

    def testFluidSystemCO2Mdot10(self):
        params.m_dot_IP = 10
        full_system.solve()

        output = full_system.gatherOutput()

        print(*testAssert(output.fluid_system_solver.pp.dP_surface, 5.9421e6, 'test_dP_surface'))
        print(*testAssert(output.fluid_system_solver.production_well.state.T_C(), 59.2802, 'test_T_prod_surface_C'))
        print(*testAssert(output.energy_results.W_net, 9.7662e4, 'test_W_net'))
        print(*testAssert(output.capital_cost_model.C_brownfield, 1.5454e7, 'test_C_brownfield_N'))
        print(*testAssert(output.capital_cost_model.C_greenfield, 3.8084e7, 'test_C_greenfield_N'))
        print(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 7.6167e-4, 'test_LCOE_brownfield'))


    def testFluidSystemCO2Mdot80(self):
        params.m_dot_IP = 80
        full_system.solve()

        output = full_system.gatherOutput()

        print(*testAssert(output.fluid_system_solver.pp.dP_surface, 5.6501e6, 'test_dP_surface'))
        print(*testAssert(output.fluid_system_solver.production_well.state.T_C(), 59.0540, 'test_T_prod_surface_C'))
        print(*testAssert(output.energy_results.W_net, 4.9940e+05, 'test_W_net'))
        print(*testAssert(output.capital_cost_model.C_brownfield, 2.9832e+07, 'test_C_brownfield_N'))
        print(*testAssert(output.capital_cost_model.C_greenfield, 5.2462e7, 'test_C_greenfield_N'))
        print(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 2.8754e-4, 'test_LCOE_brownfield'))

    def testFluidSystemCO2Mdot200(self):
        params.m_dot_IP = 200
        full_system.solve()

        output = full_system.gatherOutput()

        print(*testAssert(output.fluid_system_solver.pp.dP_surface, 3.468765e+06, 'test_dP_surface'))
        print(*testAssert(output.fluid_system_solver.production_well.state.T_C(), 47.6070, 'test_T_prod_surface_C'))
        print(*testAssert(output.energy_results.W_net, -1.3602e6, 'test_W_net'))
        print(*testAssert(output.capital_cost_model.C_brownfield, 5.26260e7, 'test_C_brownfield_N'))
        print(*testAssert(output.capital_cost_model.C_greenfield, 7.52561e7, 'test_C_greenfield_N'))
        print(np.isnan(output.capital_cost_model.LCOE_brownfield.LCOE), 'test_LCOE_brownfield')

    def testFluidSystemCO2Mdot100(self):
        params.m_dot_IP = 100
        params.depth = 2400.
        params.permeability = 1e-8 / 100.
        full_system.solve()

        output = full_system.gatherOutput()

        print(*testAssert(output.fluid_system_solver.production_well.state.T_C(), 55.3144, 'test_T_prod_surface_C'))

    def testFluidSystemCO2SolverOptMdot(self):

        full_system_solver = FullSystemSolverMinLCOEBrownfield(full_system)

        optMdot = full_system_solver.solve(time_years = 1)

        output = full_system.gatherOutput()
        print(*testAssert(optMdot, 56.9487, 'test_optMdot_solver_optMdot'))
        print(*testAssert(output.energy_results.W_net, 4.6239e5, 'test_optMdot_solver_w_net'))
        print(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 2.6207-4, 'test_optMdot_solver_LCOE_brownfield'))
