import numpy as np

from utils.readXlsxData import readCostTable
from src.capitalCostWell import CapitalCostWell

class CapitalCostExploration(object):
    """CapitalCostExploration."""
    @staticmethod
    def cO2Ideal(N, well_length, well_radius, success_rate, cost_year):
        A = cModeling(cost_year)
        B = dCCharacterizationWellsIdeal(cost_year, well_length, well_radius, success_rate)
        C = dCModelingCO2(N, cost_year)
        return A + B + C

    @staticmethod
    def cO2Baseline(N, well_length, well_radius, success_rate, cost_year):
        A = cModeling(cost_year)
        B = dCCharacterizationWellsBaseline(cost_year, well_length, well_radius, success_rate)
        C = dCModelingCO2(N, cost_year)
        return A + B + C

    @staticmethod
    def waterIdeal(well_length, well_radius, success_rate, cost_year):
        A = cModeling(cost_year)
        B = dCCharacterizationWellsIdeal(cost_year, well_length, well_radius, success_rate)
        return A + B

    @staticmethod
    def waterBaseline(well_length, well_radius, success_rate, cost_year):
        A = cModeling(cost_year)
        B = dCCharacterizationWellsBaseline(cost_year, well_length, well_radius, success_rate)
        return A + B

def aCO2AMA(N):
    # Limit N to at least 1
    return (np.maximum(N, 1) * 1000 + 1600)**2

X_PC_expl = 1.15
X_IC_expl = 1.05
g_characterization_wells = 2

def cModeling(cost_year):
    return X_IC_expl * X_PC_expl * readCostTable(cost_year, 'PPI_O&G-s') * 508000

def dCCharacterizationWellsIdeal(cost_year, well_length, well_radius, success_rate):
    c_well = CapitalCostWell.waterIdeal(well_length = well_length, well_radius = well_radius, success_rate = 1., cost_year = cost_year)
    return 0.2 * c_well * g_characterization_wells / success_rate

def dCCharacterizationWellsBaseline(cost_year, well_length, well_radius, success_rate):
    c_well = CapitalCostWell.waterBaseline(well_length = well_length, well_radius = well_radius, success_rate = 1., cost_year = cost_year)
    return 0.2 * c_well * g_characterization_wells / success_rate

def dCModelingCO2(N, cost_year):
    return X_IC_expl * X_PC_expl * readCostTable(cost_year, 'PPI_O&G-s') * 44800 * (1/1e6) * aCO2AMA(N)
