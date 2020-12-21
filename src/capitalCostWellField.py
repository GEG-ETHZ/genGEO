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

from utils.readXlsxData import readCostTable
from src.capitalCostWell import CapitalCostWell
from models.simulationParameters import SimulationParameters

class CapitalCostWellField(object):
    """CapitalCostWellField."""
    @staticmethod
    def cO2MonitoringBaseline(params = None, **kwargs):
        if params == None:
            params = SimulationParameters(**kwargs)
        return cPermitting(params.cost_year) + dCPermittingCO2(params.N_5spot, params.cost_year) + dCMonitoringCO2Baseline(params.N_5spot, params.cost_year, params.depth, params.monitoring_well_radius)

    @staticmethod
    def cO2MonitoringIdeal(params = None, **kwargs):
        if params == None:
            params = SimulationParameters(**kwargs)
        return cPermitting(params.cost_year) + dCPermittingCO2(params.N_5spot, params.cost_year) + dCMonitoringCO2Ideal(params.N_5spot, params.cost_year, params.depth, params.monitoring_well_radius)

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

def dCMonitoringCO2Ideal(N, cost_year, depth, monitoring_well_radius):
    c_surface_monitoring_CO2 = X_IC_wf * X_PC_wf * readCostTable('PPI_O&G-s', cost_year = cost_year) * 138000 * (1/1e6) * aCO2AMA(N)
    c_monitoring_Wells_CO2 = N**2 * CapitalCostWell.waterIdeal(depth = depth, well_radius = monitoring_well_radius, success_rate = 1., cost_year = cost_year)
    return c_monitoring_Wells_CO2 + c_surface_monitoring_CO2

def dCMonitoringCO2Baseline(N, cost_year, depth, monitoring_well_radius):
    c_surface_monitoring_CO2 = X_IC_wf * X_PC_wf * readCostTable('PPI_O&G-s', cost_year = cost_year) * 138000 * (1/1e6) * aCO2AMA(N)
    c_monitoring_Wells_CO2 = N**2 * CapitalCostWell.waterBaseline(depth = depth, well_radius = monitoring_well_radius, success_rate = 1., cost_year = cost_year)
    return c_monitoring_Wells_CO2 + c_surface_monitoring_CO2
