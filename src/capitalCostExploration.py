import numpy as np

from utils.readXlsxData import readCostTable
from src.capitalCostWell import CapitalCostWell

class CapitalCostExploration(object):
    """CapitalCostExploration."""
    @staticmethod
    def cO2Ideal(N, cost_year, well_length, well_diameter, success_rate_well):
        A = cModeling(cost_year)
        B = dCCharacterizationWellsIdeal(cost_year, well_length, well_diameter, success_rate_well)
        C = dCModelingCO2(N, cost_year)
        return A + B + C

    @staticmethod
    def cO2Baseline(N, cost_year, well_length, well_diameter, success_rate_well):
        A = cModeling(cost_year)
        B = dCCharacterizationWellsBaseline(cost_year, well_length, well_diameter, success_rate_well)
        C = dCModelingCO2(N, cost_year)
        return A + B + C

    @staticmethod
    def waterIdeal(cost_year, well_length, well_diameter, success_rate_well):
        A = cModeling(cost_year)
        B = dCCharacterizationWellsIdeal(cost_year, well_length, well_diameter, success_rate_well)
        return A + B

    @staticmethod
    def waterBaseline(cost_year, well_length, well_diameter, success_rate_well):
        A = cModeling(cost_year)
        B = dCCharacterizationWellsBaseline(cost_year, well_length, well_diameter, success_rate_well)
        return A + B

def aCO2AMA(N):
    # Limit N to at least 1
    return (np.maximum(N, 1) * 1000 + 1600)**2

X_PC_expl = 1.15
X_IC_expl = 1.05
g_characterization_wells = 2

def cModeling(cost_year):
    return X_IC_expl * X_PC_expl * readCostTable(cost_year, 'PPI_O&G-s') * 508000

def dCCharacterizationWellsIdeal(cost_year, well_length, well_diameter, success_rate_well):
    c_well = CapitalCostWell.waterIdeal(well_length, well_diameter, 1., cost_year)
    return 0.2 * c_well * g_characterization_wells / success_rate_well

def dCCharacterizationWellsBaseline(cost_year, well_length, well_diameter, success_rate_well):
    c_well = CapitalCostWell.waterBaseline(well_length, well_diameter, 1., cost_year)
    return 0.2 * c_well * g_characterization_wells / success_rate_well

def dCModelingCO2(N, cost_year):
    return X_IC_expl * X_PC_expl * readCostTable(cost_year, 'PPI_O&G-s') * 44800 * (1/1e6) * aCO2AMA(N)
