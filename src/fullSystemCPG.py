
from src.energyConversion import EnergyConversionCPG
from src.fullSystem import FullSystemOutput

from src.semiAnalyticalWell import SemiAnalyticalWell
from src.porousReservoir import PorousReservoir
from src.subsurfaceComponents import DownHolePump
from src.oRCCycleTboil import ORCCycleTboil

from src.fluidSystemCO2 import FluidSystemCO2
from src.capitalCostSystem import CapitalCostSystem

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

class FullSystemCPG(object):
    """FullSystemCPG."""

    def __init__(self, params, fluid_system_solver, capital_cost_model):
        self.params = params
        self.fluid_system_solver = fluid_system_solver
        self.capital_cost_model = capital_cost_model

    def solve(self):

        results = FullSystemOutput()

        results.fluid_system_solver = self.fluid_system_solver.solve()
        results.energy_results = EnergyConversionCPG.compute(self.params, results.fluid_system_solver)
        self.capital_cost_model.energy_results = results.energy_results
        self.capital_cost_model.fluid_system = results.fluid_system_solver
        results.capital_cost_model = self.capital_cost_model.solve()

        return results

    #static method below instantiates a default CPG system
    @staticmethod
    def getDefaultCPGSystem(params = None, **kwargs):
        if params == None:
            params = SimulationParameters(**kwargs)

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

        return FullSystemCPG(params, fluid_system, capital_cost_system)
