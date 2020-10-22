
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

from src.totalSystemBase import TotalSystemBase

from utils.simulationParameters import SimulationParameters

class TotalSystemCO2(TotalSystemBase):
    """TotalSystemCO2."""

    def __init__(self, **kwargs):

        # define global methods to be used in this tests
        gsp = SimulationParameters(working_fluid = 'co2', **kwargs)

        fluid_system = FluidSystemCO2(T_ambient_C = gsp.T_ambient_C,
                                dT_approach = gsp.dT_approach,
                                eta_pump = gsp.eta_pump_co2,
                                eta_turbine = gsp.eta_turbine_co2,
                                coolingMode = gsp.cooling_mode)

        fluid_system.injection_well = SemiAnalyticalWell(params = gsp,
                                            dz_total = -gsp.depth,
                                            wellRadius = gsp.well_radius,
                                            wellMultiplier = gsp.well_multiplier,
                                            fluid = gsp.working_fluid,
                                            epsilon = gsp.epsilon,
                                            dT_dz = gsp.dT_dz,
                                            T_e_initial = gsp.T_e_initial)

        fluid_system.reservoir = PorousReservoir(params = gsp)

        fluid_system.production_well = SemiAnalyticalWell(params = gsp,
                                            dz_total = gsp.depth,
                                            wellRadius = gsp.well_radius,
                                            wellMultiplier = gsp.well_multiplier,
                                            fluid = gsp.working_fluid,
                                            epsilon = gsp.epsilon,
                                            dT_dz = gsp.dT_dz,
                                            T_e_initial = gsp.T_e_initial + gsp.dT_dz * gsp.depth)


        capital_cost_system = CapitalCostSystem(gsp.cost_year)
        capital_cost_system.CapitalCost_SurfacePlant = CapitalCostSurfacePlantCPG(gsp.cost_year)
        capital_cost_system.CapitalCost_SurfacePipe = CapitalCostSurfacePipes(N = 0)
        capital_cost_system.CapitalCost_Production_Well = CapitalCostWell.cO2Baseline(well_length = gsp.depth,
                                                                            well_diameter = gsp.well_radius * 2,
                                                                            success_rate = gsp.well_success_rate,
                                                                            cost_year = gsp.cost_year)
        capital_cost_system.CapitalCost_Injection_Well = CapitalCostWell.cO2Baseline(well_length = gsp.depth,
                                                                            well_diameter = gsp.well_radius * 2,
                                                                            success_rate = gsp.well_success_rate,
                                                                            cost_year = gsp.cost_year)
        capital_cost_system.CapitalCost_Wellfield = CapitalCostWellField.cO2MonitoringBaseline(N=1,
                                                                            monitoring_well_length = gsp.depth,
                                                                            monitoring_well_diameter = gsp.monitoring_well_diameter,
                                                                            cost_year=gsp.cost_year)
        capital_cost_system.CapitalCost_Exploration = CapitalCostExploration.cO2Baseline(N=1,
                                                                            well_length = gsp.depth,
                                                                            well_diameter = gsp.well_radius * 2,
                                                                            success_rate = gsp.well_success_rate,
                                                                            cost_year = gsp.cost_year)
        capital_cost_system.lcoe_model = LCOESimple(F_OM = gsp.F_OM,
                                                    discountRate = gsp.discount_rate,
                                                    Lifetime = gsp.lifetime,
                                                    CapacityFactor = gsp.capacity_factor)

        self.full_system = FullSystemCPG(fluid_system, capital_cost_system)
