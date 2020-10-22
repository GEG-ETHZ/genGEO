
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

from src.totalSystemBase import TotalSystemBase

from utils.simulationParameters import SimulationParameters

class TotalSystemWater(TotalSystemBase):
    """TotalSystemWater."""

    def __init__(self, **kwargs):

        # define global methods to be used in this tests
        gsp = SimulationParameters(working_fluid = 'water', **kwargs)

        fluid_system = FluidSystemWater()
        fluid_system.injection_well = SemiAnalyticalWell(params = gsp,
                                            dz_total = -1 * gsp.depth,
                                            wellRadius = gsp.well_radius,
                                            wellMultiplier = gsp.well_multiplier,
                                            fluid = gsp.working_fluid,
                                            epsilon = gsp.epsilon,
                                            dT_dz = gsp.dT_dz,
                                            T_e_initial = gsp.T_e_initial)

        fluid_system.reservoir = PorousReservoir(params = gsp)

        fluid_system.production_well1 = SemiAnalyticalWell(params = gsp,
                                            dz_total = gsp.depth - gsp.pump_depth,
                                            wellRadius = gsp.well_radius,
                                            wellMultiplier = gsp.well_multiplier,
                                            fluid = gsp.working_fluid,
                                            epsilon = gsp.epsilon,
                                            dT_dz = gsp.dT_dz,
                                            T_e_initial = gsp.T_e_initial + gsp.dT_dz * gsp.depth)

        prod_well2 = SemiAnalyticalWell(params = gsp,
                                            dz_total = gsp.pump_depth,
                                            wellRadius = gsp.well_radius,
                                            wellMultiplier = gsp.well_multiplier,
                                            fluid = gsp.working_fluid,
                                            epsilon = gsp.epsilon,
                                            dT_dz = gsp.dT_dz,
                                            T_e_initial = gsp.T_e_initial + gsp.dT_dz * gsp.pump_depth)

        fluid_system.pump = DownHolePump(well = prod_well2,
                            pump_depth = gsp.pump_depth,
                            max_pump_dP = gsp.max_pump_dP,
                            eta_pump = gsp.eta_pump)

        fluid_system.pp = ORCCycleTboil(T_ambient_C = gsp.T_ambient_C,
                                dT_approach = gsp.dT_approach,
                                dT_pinch = gsp.dT_pinch,
                                eta_pump = gsp.eta_pump_orc,
                                eta_turbine = gsp.eta_turbine_orc,
                                coolingMode = gsp.cooling_mode,
                                orcFluid = gsp.orc_fluid)

        capital_cost_system = CapitalCostSystem(gsp.cost_year)
        capital_cost_system.CapitalCost_SurfacePlant = CapitalCostSurfacePlantORC(gsp.cost_year)
        capital_cost_system.CapitalCost_SurfacePipe = CapitalCostSurfacePipes(N = 0)
        capital_cost_system.CapitalCost_Production_Well = CapitalCostWell.waterBaseline(well_length = gsp.depth,
                                                                            well_diameter = gsp.well_radius * 2,
                                                                            success_rate = gsp.well_success_rate,
                                                                            cost_year = gsp.cost_year)
        capital_cost_system.CapitalCost_Injection_Well = CapitalCostWell.waterBaseline(well_length = gsp.depth,
                                                                            well_diameter = gsp.well_radius * 2,
                                                                            success_rate = gsp.well_success_rate,
                                                                            cost_year = gsp.cost_year)
        capital_cost_system.CapitalCost_Wellfield = CapitalCostWellField.water(gsp.cost_year)
        capital_cost_system.CapitalCost_Exploration = CapitalCostExploration.waterBaseline(well_length = gsp.depth,
                                                                            well_diameter = gsp.well_radius * 2,
                                                                            success_rate = gsp.well_success_rate,
                                                                            cost_year = gsp.cost_year)
        capital_cost_system.lcoe_model = LCOESimple(F_OM = gsp.F_OM,
                                                    discountRate = gsp.discount_rate,
                                                    Lifetime = gsp.lifetime,
                                                    CapacityFactor = gsp.capacity_factor)

        solver = FluidSystemWaterSolver(fluid_system)

        self.full_system = FullSystemORC(solver, capital_cost_system)
