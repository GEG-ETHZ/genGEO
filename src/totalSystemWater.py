
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
from src.capitalCostSurfacePlantORC import CapitalCostSurfacePlantORC
from src.capitalCostWellStimulation import CapitalCostWellStimulation

from src.totalSystemBase import TotalSystemBase

from utils.simulationParameters import SimulationParameters

class TotalSystemWater(TotalSystemBase):
    """TotalSystemWater."""

    def __init__(self, **kwargs):

        # define global methods to be used in this tests
        params = SimulationParameters(working_fluid = 'water', **kwargs)

        fluid_system = FluidSystemWater(params = params)
        fluid_system.injection_well = SemiAnalyticalWell(params = params,
                                            dz_total = -1 * params.depth,
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

        solver = FluidSystemWaterSolver(fluid_system)

        self.full_system = FullSystemORC(params, solver, capital_cost_system)
