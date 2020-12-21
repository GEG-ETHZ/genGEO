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
import numpy as np

from src.capitalCostWell import CapitalCostWell

from utils.readXlsxData import readCostTable
from models.simulationParameters import SimulationParameters

class CapitalCostExploration(object):
    """CapitalCostExploration."""
    @staticmethod
    def cO2Ideal(params = None, **kwargs):
        if params == None:
            params = SimulationParameters(**kwargs)
        A = cModeling(params.cost_year)
        B = dCCharacterizationWellsIdeal(params.cost_year, params.depth, params.well_radius, params.success_rate)
        C = dCModelingCO2(params.N_5spot, params.cost_year)
        return A + B + C

    @staticmethod
    def cO2Baseline(params = None, **kwargs):
        if params == None:
            params = SimulationParameters(**kwargs)
        A = cModeling(params.cost_year)
        B = dCCharacterizationWellsBaseline(params.cost_year, params.depth, params.well_radius, params.success_rate)
        C = dCModelingCO2(params.N_5spot, params.cost_year)
        return A + B + C

    @staticmethod
    def waterIdeal(params = None, **kwargs):
        if params == None:
            params = SimulationParameters(**kwargs)
        A = cModeling(params.cost_year)
        B = dCCharacterizationWellsIdeal(params.cost_year, params.depth, params.well_radius, params.success_rate)
        return A + B

    @staticmethod
    def waterBaseline(params = None, **kwargs):
        if params == None:
            params = SimulationParameters(**kwargs)
        A = cModeling(params.cost_year)
        B = dCCharacterizationWellsBaseline(params.cost_year, params.depth, params.well_radius, params.success_rate)
        return A + B

def aCO2AMA(N):
    # Limit N to at least 1
    return (np.maximum(N, 1) * 1000 + 1600)**2

X_PC_expl = 1.15
X_IC_expl = 1.05
g_characterization_wells = 2

def cModeling(cost_year):
    return X_IC_expl * X_PC_expl * readCostTable('PPI_O&G-s', cost_year = cost_year) * 508000

def dCCharacterizationWellsIdeal(cost_year, depth, well_radius, success_rate):
    c_well = CapitalCostWell.waterIdeal(depth = depth, well_radius = well_radius, success_rate = 1., cost_year = cost_year)
    return 0.2 * c_well * g_characterization_wells / success_rate

def dCCharacterizationWellsBaseline(cost_year, depth, well_radius, success_rate):
    c_well = CapitalCostWell.waterBaseline(depth = depth, well_radius = well_radius, success_rate = 1., cost_year = cost_year)
    return 0.2 * c_well * g_characterization_wells / success_rate

def dCModelingCO2(N, cost_year):
    return X_IC_expl * X_PC_expl * readCostTable('PPI_O&G-s', cost_year = cost_year) * 44800 * (1/1e6) * aCO2AMA(N)
