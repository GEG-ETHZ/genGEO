import numpy as np

from utils.readXlsxData import readCostTable
from src.capitalCostWell import CapitalCostWell


class CapitalCostWellField(object):
    """CapitalCostWellField."""
    @staticmethod
    def cO2MonitoringBaseline(params = None, **kwargs):
        if params == None:
            params = SimulationParameters(**kwargs)
        return cPermitting(params.cost_year) + dCPermittingCO2(params.well_cost_N, params.cost_year) + dCMonitoringCO2Baseline(params.well_cost_N, params.cost_year, params.monitoring_well_length, params.monitoring_well_radius)

    @staticmethod
    def cO2MonitoringIdeal(params = None, **kwargs):
        if params == None:
            params = SimulationParameters(**kwargs)
        return cPermitting(params.cost_year) + dCPermittingCO2(params.well_cost_N, params.cost_year) + dCMonitoringCO2Ideal(params.well_cost_N, params.cost_year, params.monitoring_well_length, params.monitoring_well_radius)

    @staticmethod
    def cO2(params = None, **kwargs):
        if params == None:
            params = SimulationParameters(**kwargs)
        return cPermitting(params.cost_year)

    @staticmethod
    def water(params = None, **kwargs):
        if params == None:
            params = SimulationParameters(**kwargs)
        return cPermitting(params.cost_year)

X_IC_wf = 1.05
X_PC_wf = 1.15

def aCO2AMA(N):
    # Limit N to at least 1
    return (np.maximum(N, 1) * 1000 + 1600)**2

def cPermitting(cost_year):
    return X_IC_wf * X_PC_wf * readCostTable('PPI_Permit', cost_year = cost_year) * 665700.

def dCPermittingCO2(N, cost_year):
    return X_IC_wf * X_PC_wf * readCostTable('PPI_Permit', cost_year = cost_year) * 45000 * (1/1e6) * aCO2AMA(N)

def dCMonitoringCO2Ideal(N, cost_year, monitoring_well_length, monitoring_well_radius):
    c_surface_monitoring_CO2 = X_IC_wf * X_PC_wf * readCostTable('PPI_O&G-s', cost_year = cost_year) * 138000 * (1/1e6) * aCO2AMA(N)
    c_monitoring_Wells_CO2 = N**2 * CapitalCostWell.waterIdeal(well_length = monitoring_well_length, well_radius = monitoring_well_radius, success_rate = 1., cost_year = cost_year)
    return c_monitoring_Wells_CO2 + c_surface_monitoring_CO2

def dCMonitoringCO2Baseline(N, cost_year, monitoring_well_length, monitoring_well_radius):
    c_surface_monitoring_CO2 = X_IC_wf * X_PC_wf * readCostTable('PPI_O&G-s', cost_year = cost_year) * 138000 * (1/1e6) * aCO2AMA(N)
    c_monitoring_Wells_CO2 = N**2 * CapitalCostWell.waterBaseline(well_length = monitoring_well_length, well_radius = monitoring_well_radius, success_rate = 1., cost_year = cost_year)
    return c_monitoring_Wells_CO2 + c_surface_monitoring_CO2
