
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
from src.capitalCostSurfacePlantCPG import CapitalCostSurfacePlantCPG
from src.capitalCostWellStimulation import CapitalCostWellStimulation

from src.totalSystemBase import TotalSystemBase

from utils.simulationParameters import SimulationParameters

class TotalSystemCO2(TotalSystemBase):
    """TotalSystemCO2."""

    def __init__(self, **kwargs):

        params = SimulationParameters(working_fluid = 'co2', **kwargs)

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

        self.full_system = FullSystemCPG(params, fluid_system, capital_cost_system)
