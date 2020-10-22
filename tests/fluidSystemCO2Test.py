import unittest
import numpy as np

from src.semiAnalyticalWell import SemiAnalyticalWell
from src.porousReservoir import PorousReservoir
from src.subsurfaceComponents import DownHolePump
from src.oRCCycleTboil import ORCCycleTboil

from src.fluidSystemCO2 import FluidSystemCO2
from src.capitalCostSystem import CapitalCostSystem

from src.fullSystem import FullSystemCPG
from src.lCOESimple import LCOESimple
from src.capitalCostSurfacePipes import CapitalCostSurfacePipes
from src.capitalCostWell import CapitalCostWell
from src.capitalCostWellField import CapitalCostWellField
from src.capitalCostExploration import CapitalCostExploration
from src.capitalCostSurfacePlantCPG import CapitalCostSurfacePlantCPG
from src.fullSystemSolver import FullSystemSolverMinLCOEBrownfield

from utils.simulationParameters import SimulationParameters
from utils.fluidStateFromPT import FluidStateFromPT

from tests.testAssertion import testAssert


# load simulaiton parameters
params = SimulationParameters(working_fluid = 'co2')
params.wellMultiplier = 4.

fluid_system = FluidSystemCO2(params = params)

fluid_system.injection_well = SemiAnalyticalWell(params = params,
                                    dz_total = -params.depth,
                                    T_e_initial = params.T_ambient_C)

fluid_system.reservoir = PorousReservoir(params = params)

fluid_system.production_well = SemiAnalyticalWell(params = params,
                                    dz_total = params.depth,
                                    T_e_initial = params.T_ambient_C + params.dT_dz * params.depth)


capital_cost_system = CapitalCostSystem(params)
capital_cost_system.CapitalCost_SurfacePlant = CapitalCostSurfacePlantCPG(params)
capital_cost_system.CapitalCost_SurfacePipe = CapitalCostSurfacePipes(N = 0)
capital_cost_system.CapitalCost_Production_Well = CapitalCostWell.cO2Baseline(params)
capital_cost_system.CapitalCost_Injection_Well = CapitalCostWell.cO2Baseline(params)
capital_cost_system.CapitalCost_Wellfield = CapitalCostWellField.cO2MonitoringBaseline(N=1,
                                                                    monitoring_well_length=2500,
                                                                    monitoring_well_diameter=0.108,
                                                                    cost_year=2019)
capital_cost_system.CapitalCost_Exploration = CapitalCostExploration.cO2Baseline(N=1,
                                                                    well_length = 2500,
                                                                    well_diameter = 0.205,
                                                                    success_rate = 0.95,
                                                                    cost_year = 2019)
capital_cost_system.lcoe_model = LCOESimple(F_OM = 0.045,
                                            discountRate = 0.096,
                                            Lifetime = 25,
                                            CapacityFactor = 0.9)

full_system = FullSystemCPG(params, fluid_system, capital_cost_system)

class FluidSystemCO2Test(unittest.TestCase):

    def testFluidSystemCO2Mdot10(self):
        full_system.solve(m_dot = 10, time_years = 1)

        output = full_system.gatherOutput()

        self.assertTrue(*testAssert(output.fluid_system_solver.pp.dP_surface, 6.3943e6, 'test_dP_surface'))
        self.assertTrue(*testAssert(output.fluid_system_solver.production_well.end_T_C(), 60.4516, 'test_T_prod_surface_C'))
        self.assertTrue(*testAssert(output.energy_results.W_net, 1.0283e5, 'test_W_net'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_brownfield, 1.5452e7, 'test_C_brownfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_greenfield, 3.8082e7, 'test_C_greenfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 7.2330e-4, 'test_LCOE_brownfield'))


    def testFluidSystemCO2Mdot80(self):
        full_system.solve(m_dot = 80, time_years = 1)

        output = full_system.gatherOutput()

        self.assertTrue(*testAssert(output.fluid_system_solver.pp.dP_surface, 5.65949586e6, 'test_dP_surface'))
        self.assertTrue(*testAssert(output.fluid_system_solver.production_well.end_T_C(), 59.0826, 'test_T_prod_surface_C'))
        self.assertTrue(*testAssert(output.energy_results.W_net, 5.02286e+05, 'test_W_net'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_brownfield, 2.98408178e+07, 'test_C_brownfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_greenfield, 5.2470e7, 'test_C_greenfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 2.8596677e-04, 'test_LCOE_brownfield'))

    def testFluidSystemCO2Mdot200(self):
        full_system.solve(m_dot = 200, time_years = 1)

        output = full_system.gatherOutput()

        self.assertTrue(*testAssert(output.fluid_system_solver.pp.dP_surface, 3.468765e+06, 'test_dP_surface'))
        self.assertTrue(*testAssert(output.fluid_system_solver.production_well.end_T_C(), 47.5801, 'test_T_prod_surface_C'))
        self.assertTrue(*testAssert(output.energy_results.W_net, -1.3602e6, 'test_W_net'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_brownfield, 5.26260e7, 'test_C_brownfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_greenfield, 7.52561e7, 'test_C_greenfield_N'))
        self.assertTrue(np.isnan(output.capital_cost_model.LCOE_brownfield.LCOE), 'test_LCOE_brownfield')

    def testFluidSystemCO2SolverOptMdot(self):

        full_system_solver = FullSystemSolverMinLCOEBrownfield(full_system)

        optMdot = full_system_solver.solve(time_years = 1)

        output = full_system.gatherOutput()
        print(*testAssert(optMdot, 56.9487, 'test_optMdot_solver_optMdot'))
        print(*testAssert(output.energy_results.W_net, 4.6239e5, 'test_optMdot_solver_w_net'))
        print(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 2.6207-4, 'test_optMdot_solver_LCOE_brownfield'))
