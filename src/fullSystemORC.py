# Licensed under LGPL 2.1, please see LICENSE for details
# https://www.gnu.org/licenses/lgpl-2.1.html
#
# The work on this project has been performed at the GEG Group at ETH Zurich:
# --> https://geg.ethz.ch
#
# The initial version of this file has been implemented by:
#
#     Philipp Schaedle (https://github.com/philippschaedle)
#     Benjamin M. Adams
#
# Further changes are done by:
#

############################
from src.energyConversion import EnergyConversionORC
from src.fullSystem import FullSystemOutput

from src.semiAnalyticalWell import SemiAnalyticalWell
from src.porousReservoir import PorousReservoir
from src.downholePump import DownHolePump
from src.oRCCycleTboil import ORCCycleTboil

from src.fluidSystemWater import FluidSystemWater
from src.fluidSystemWaterSolver import FluidSystemWaterSolver
from src.capitalCostSystem import CapitalCostSystem

from src.lCOESimple import LCOESimple
from src.capitalCostSurfacePipes import CapitalCostSurfacePipes
from src.capitalCostWell import CapitalCostWell
from src.capitalCostWellField import CapitalCostWellField
from src.capitalCostExploration import CapitalCostExploration
from src.capitalCostWellStimulation import CapitalCostWellStimulation
from src.capitalCostSurfacePlantORC import CapitalCostSurfacePlantORC

from models.simulationParameters import SimulationParameters

class FullSystemORC(object):
    """FullSystemORC."""

    def __init__(self, params, fluid_system_solver, capital_cost_model):
        self.params = params
        self.fluid_system_solver = fluid_system_solver
        self.capital_cost_model = capital_cost_model

    def solve(self):

        results = FullSystemOutput()

        results.fluid_system_solver = self.fluid_system_solver.solve()
        results.energy_results = EnergyConversionORC.compute(self.params, results.fluid_system_solver)
        self.capital_cost_model.energy_results = results.energy_results
        self.capital_cost_model.fluid_system = results.fluid_system_solver
        results.capital_cost_model = self.capital_cost_model.solve()

        return results

    def getDefaultWaterSystem(params = None, **kwargs):
        if params == None:
            params = SimulationParameters(**kwargs)

        fluid_system = FluidSystemWater(params = params)
        injWell_m_dot_multiplier = params.wellFieldType.getInjWellMdotMultiplier()
        fluid_system.injection_well = SemiAnalyticalWell(params = params,
                                            dz_total = -params.depth,
                                            T_e_initial = params.T_ambient_C,
                                            m_dot_multiplier = injWell_m_dot_multiplier)
        fluid_system.reservoir = PorousReservoir(params = params)
        prodWell_m_dot_multiplier = params.wellFieldType.getProdWellMdotMultiplier()
        fluid_system.production_well1 = SemiAnalyticalWell(params = params,
                                            dz_total = params.depth - params.pump_depth,
                                            T_e_initial = params.T_ambient_C + params.dT_dz * params.depth,
                                            m_dot_multiplier = prodWell_m_dot_multiplier)
        prod_well2 = SemiAnalyticalWell(params = params,
                                        dz_total = params.pump_depth,
                                        T_e_initial = params.T_ambient_C + params.dT_dz * params.pump_depth,
                                        m_dot_multiplier = prodWell_m_dot_multiplier)
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

        return FullSystemORC(params, solver, capital_cost_system)

