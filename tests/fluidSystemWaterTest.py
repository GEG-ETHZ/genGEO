import unittest

from src.semiAnalyticalWell import SemiAnalyticalWell
from src.porousReservoir import PorousReservoir
from src.subsurfaceComponents import DownHolePump
from src.oRCCycleTboil import ORCCycleTboil

from src.fluidSystemWater import FluidSystemWater
from src.fluidSystemWaterSolver import FluidSystemWaterSolver
from src.capitalCostSystem import CapitalCostSystem

from src.fullSystemORC import FullSystemORC
from src.lCOESimple import LCOESimple
from src.capitalCostSurfacePipes import CapitalCostSurfacePipes
from src.capitalCostWell import CapitalCostWell
from src.capitalCostWellField import CapitalCostWellField
from src.capitalCostExploration import CapitalCostExploration
from src.capitalCostWellStimulation import CapitalCostWellStimulation
from src.capitalCostSurfacePlantORC import CapitalCostSurfacePlantORC
from src.fullSystemSolver import FullSystemSolverMinLCOEBrownfield

from utils.simulationParameters import SimulationParameters
from utils.fluidStateFromPT import FluidStateFromPT

from tests.testAssertion import testAssert


# load simulaiton parameters
params = SimulationParameters(working_fluid = 'water')
params.well_multiplier = 4.
params.orc_fluid = 'R245fa'
params.capacity_factor = 0.9


fluid_system = FluidSystemWater(params = params)
fluid_system.injection_well = SemiAnalyticalWell(params = params,
                                    dz_total = -params.depth,
                                    T_e_initial = params.T_ambient_C)
fluid_system.reservoir = PorousReservoir(params = params)
fluid_system.production_well1 = SemiAnalyticalWell(params = params,
                                    dz_total = params.depth - params.pump_depth,
                                    T_e_initial = params.T_ambient_C + params.dT_dz * params.depth)
prod_well2 = SemiAnalyticalWell(params = params,
                                dz_total = params.pump_depth,
                                T_e_initial = params.T_ambient_C + params.dT_dz * params.pump_depth)
fluid_system.pump = DownHolePump(well = prod_well2,
                                params = params)
fluid_system.pp = ORCCycleTboil(params = params)

capital_cost_system = CapitalCostSystem()
capital_cost_system.CapitalCost_SurfacePlant = CapitalCostSurfacePlantORC(params = params)
capital_cost_system.CapitalCost_SurfacePipe = CapitalCostSurfacePipes.cost(params = params)
capital_cost_system.CapitalCost_Production_Well = CapitalCostWell.waterBaseline(params = params)
capital_cost_system.CapitalCost_Injection_Well = CapitalCostWell.waterBaseline(params = params)
capital_cost_system.CapitalCost_Wellfield = CapitalCostWellField.water(params = params)
capital_cost_system.CapitalCost_Exploration = CapitalCostExploration.waterBaseline(params = params)
capital_cost_system.CapitalCost_Stimulation = CapitalCostWellStimulation.cost()
capital_cost_system.lcoe_model = LCOESimple(params = params)

class FluidSystemWaterTest(unittest.TestCase):

    def testFluidSystemWaterSolverMdot1(self):
        solver = FluidSystemWaterSolver(fluid_system)

        params.m_dot_IP = 1
        full_system = FullSystemORC(params, solver, capital_cost_system)
        output = full_system.solve()

        print(*testAssert(output.fluid_system_solver.pump.well.state.P_Pa(), 7.84791e5, 'test_subsurface_solver1_pressure'))
        print(*testAssert(output.fluid_system_solver.pump.well.state.T_C(), 81.2595, 'test_subsurface_solver1_temp'))
        print(*testAssert(output.energy_results.W_net, 5.2775e3, 'test_subsurface_solver1_w_net'))
        print(*testAssert(output.capital_cost_model.C_brownfield, 9.1965e6, 'test_solver1_C_brownfield_N'))
        print(*testAssert(output.capital_cost_model.C_greenfield, 2.2308e7, 'test_solver1_C_greenfield_N'))
        print(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 0.0083879, 'test_solver1_LCOE_brownfield'))

    def testFluidSystemWaterSolverMdot40(self):
        solver = FluidSystemWaterSolver(fluid_system)

        params.m_dot_IP = 40
        full_system = FullSystemORC(params, solver, capital_cost_system)
        output = full_system.solve()

        print(*testAssert(output.fluid_system_solver.pump.well.state.P_Pa(), 7.9610e6, 'test_subsurface_solver2_pressure'))
        print(*testAssert(output.fluid_system_solver.pump.well.state.T_C(), 100.3125, 'test_subsurface_solver2_temp'))
        print(*testAssert(output.energy_results.W_net, 1.4286e+05, 'test_subsurface_solver2_w_net'))
        print(*testAssert(output.capital_cost_model.C_brownfield, 2.3287e7, 'test_solver2_C_brownfield_N'))
        print(*testAssert(output.capital_cost_model.C_greenfield, 3.6399062e7, 'test_solver2_C_greenfield_N'))
        print(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 7.8465e-4, 'test_solver2_LCOE_brownfield'))

    def testFluidSystemWaterSolverOptMdot(self):
        solver = FluidSystemWaterSolver(fluid_system)

        full_system = FullSystemORC(params, solver, capital_cost_system)

        full_system_solver = FullSystemSolverMinLCOEBrownfield(full_system)

        output = full_system_solver.solve()

        print(*testAssert(output.optMdot, 21.9823, 'test_optMdot_solver_optMdot'))
        self.assertTrue(*testAssert(output.energy_results.W_net, 1.5212e5, 'test_optMdot_solver_w_net'))
        self.assertTrue(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 5.4633e-4, 'test_optMdot_solver_brownfield'))
