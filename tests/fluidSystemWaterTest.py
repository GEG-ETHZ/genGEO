import unittest

from src.semiAnalyticalWell import SemiAnalyticalWell
from src.porousReservoir import PorousReservoir
from src.subsurfaceComponents import DownHolePump
from src.oRCCycleTboil import ORCCycleTboil

from src.fluidSystemWater import FluidSystemWater
from src.fluidSystemWaterSolver import FluidSystemWaterSolver
from src.capitalCostSystem import CapitalCostSystem

from src.fullSystem import FullSystemORC
from src.lCOESimple import LCOESimple
from src.capitalCostSurfacePipes import CapitalCostSurfacePipes
from src.capitalCostWell import CapitalCostWell
from src.capitalCostWellField import CapitalCostWellField
from src.capitalCostExploration import CapitalCostExploration
from src.capitalCostSurfacePlantORC import CapitalCostSurfacePlantORC
from src.fullSystemSolver import FullSystemSolverMinLCOEBrownfield

from utils.simulationParameters import SimulationParameters
from utils.fluidStateFromPT import FluidStateFromPT

from tests.testAssertion import testAssert


# define global methods to be used in this tests
gsp = SimulationParameters()
# alternative simulation properties for porous reservoir
gsp2 = SimulationParameters()
gsp2.k_rock = 2.1        #W/m/C
gsp2.rho_rock = 2300     #kg/m^3
gsp2.c_rock = 920.       #J/kg/C
gsp2.working_fluid = 'water'

fluid_system = FluidSystemWater()
fluid_system.injection_well = SemiAnalyticalWell(params = gsp,
                                    dz_total = -2500.,
                                    wellRadius = 0.205,
                                    wellMultiplier = 4.,
                                    fluid = 'Water',
                                    epsilon = 55 * 1e-6,
                                    dT_dz = 0.035,
                                    T_e_initial = 15.)

fluid_system.reservoir = PorousReservoir(params = gsp2)

fluid_system.production_well1 = SemiAnalyticalWell(params = gsp,
                                    dz_total = 2000.,
                                    wellRadius = 0.205,
                                    wellMultiplier = 4.,
                                    fluid = 'Water',
                                    epsilon = 55 * 1e-6,
                                    dT_dz = 0.035,
                                    T_e_initial = 15. + 0.035 * 2500.)

prod_well2 = SemiAnalyticalWell(params = gsp,
                                    dz_total = 500.,
                                    wellRadius = 0.205,
                                    wellMultiplier = 4.,
                                    fluid = 'Water',
                                    epsilon = 55 * 1e-6,
                                    dT_dz = 0.035,
                                    T_e_initial = 15. + 0.035 * 500.)

fluid_system.pump = DownHolePump(well = prod_well2,
                    pump_depth = 500.,
                    max_pump_dP = 10.e6,
                    eta_pump = 0.75)

fluid_system.pp = ORCCycleTboil(T_ambient_C = 15.,
                        dT_approach = 7.,
                        dT_pinch = 5.,
                        eta_pump = 0.9,
                        eta_turbine = 0.8,
                        coolingMode = 'Wet',
                        orcFluid = 'R245fa')

capital_cost_system = CapitalCostSystem(2019)
capital_cost_system.CapitalCost_SurfacePlant = CapitalCostSurfacePlantORC(2019)
capital_cost_system.CapitalCost_SurfacePipe = CapitalCostSurfacePipes(N = 0)
capital_cost_system.CapitalCost_Production_Well = CapitalCostWell.waterBaseline(well_length = 2500,
                                                                    well_diameter = 0.205 * 2,
                                                                    success_rate = 0.95,
                                                                    cost_year = 2019)
capital_cost_system.CapitalCost_Injection_Well = CapitalCostWell.waterBaseline(well_length = 2500,
                                                                    well_diameter = 0.205 * 2,
                                                                    success_rate = 0.95,
                                                                    cost_year = 2019)
capital_cost_system.CapitalCost_Wellfield = CapitalCostWellField.water(2019)
capital_cost_system.CapitalCost_Exploration = CapitalCostExploration.waterBaseline(well_length = 2500,
                                                                    well_diameter = 0.205 * 2,
                                                                    success_rate = 0.95,
                                                                    cost_year = 2019)
capital_cost_system.lcoe_model = LCOESimple(F_OM = 0.045,
                                            discountRate = 0.096,
                                            Lifetime = 25,
                                            CapacityFactor = 0.9)

initialState = FluidStateFromPT(1.e6, 60., fluid_system.fluid)

class FluidSystemWaterTest(unittest.TestCase):

    def testFluidSystemWaterMdot1(self):
        fluid_system_results = fluid_system.solve(initial_state = initialState,
                        m_dot = 1,
                        time_years = 1)

        self.assertTrue(*testAssert(fluid_system_results.P_Pa(), 1.352117133629285e+06, 'test_subsurface1_pressure'))
        self.assertTrue(*testAssert(fluid_system_results.T_C(), 54.038702389379260, 'test_subsurface1_temp'))

    def testFluidSystemWaterSolverMdot1(self):
        solver = FluidSystemWaterSolver(fluid_system)

        full_system = FullSystemORC(solver, capital_cost_system)
        full_system.solve(m_dot = 1, time_years = 1)

        output = full_system.gatherOutput()

        self.assertTrue(*testAssert(output.fluid_system_solver.production_well2.end_P_Pa(), 1.400500841642725e+06, 'test_subsurface_solver1_pressure'))
        self.assertTrue(*testAssert(output.fluid_system_solver.production_well2.end_T_C(), 81.255435963675800, 'test_subsurface_solver1_temp'))
        self.assertTrue(*testAssert(output.energy_results.W_net, 5.3517e3, 'test_subsurface_solver1_w_net'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_brownfield, 1.1467e7, 'test_solver1_C_brownfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_greenfield, 2.4579e7, 'test_solver1_C_greenfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 0.010314, 'test_solver1_LCOE_brownfield'))

    def testFluidSystemWaterMdot40(self):
        fluid_system_results = fluid_system.solve(initial_state = initialState,
                        m_dot = 40,
                        time_years = 1)

        self.assertTrue(*testAssert(fluid_system_results.P_Pa(), 8.149344891612666e+06, 'test_subsurface2_pressure'))
        self.assertTrue(*testAssert(fluid_system_results.T_C(), 58.835063067362940, 'test_subsurface2_temp'))

    def testFluidSystemWaterSolverMdot40(self):
        solver = FluidSystemWaterSolver(fluid_system)

        full_system = FullSystemORC(solver, capital_cost_system)
        full_system.solve(m_dot = 40, time_years = 1)

        output = full_system.gatherOutput()

        self.assertTrue(*testAssert(output.fluid_system_solver.production_well2.end_P_Pa(), 8.1493e+06, 'test_subsurface_solver2_pressure'))
        self.assertTrue(*testAssert(output.fluid_system_solver.production_well2.end_T_C(), 100.3127, 'test_subsurface_solver2_temp'))
        self.assertTrue(*testAssert(output.energy_results.W_net, 1.53849e5, 'test_subsurface_solver2_w_net'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_brownfield, 2.6733e7, 'test_solver2_C_brownfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_greenfield, 3.9845e7, 'test_solver2_C_greenfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 8.36432e-04, 'test_solver2_LCOE_brownfield'))

    def testFluidSystemWaterSolverOptMdot(self):
        solver = FluidSystemWaterSolver(fluid_system)

        full_system = FullSystemORC(solver, capital_cost_system)

        full_system_solver = FullSystemSolverMinLCOEBrownfield(full_system)

        optMdot = full_system_solver.solve(time_years = 1)

        output = full_system.gatherOutput()
        print(*testAssert(optMdot, 22.85, 'test_optMdot_solver_optMdot'))
        self.assertTrue(*testAssert(output.energy_results.W_net, 1.6292e5, 'test_optMdot_solver_w_net'))
        self.assertTrue(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 6.0634e-4, 'test_optMdot_solver_brownfield'))
