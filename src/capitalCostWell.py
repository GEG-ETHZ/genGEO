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

from utils.readXlsxData import readCostTable
from models.simulationParameters import SimulationParameters


class CapitalCostWell(object):
    """CapitalCostWell."""
    @staticmethod
    def cO2Baseline(params = None, **kwargs):
        if params == None:
            params = SimulationParameters(**kwargs)
        return wellCO2(params.depth, params.well_radius, baseLine(params.depth), params.success_rate, params.cost_year)

    @staticmethod
    def cO2Ideal(params = None, **kwargs):
        if params == None:
            params = SimulationParameters(**kwargs)
        return wellCO2(params.depth, params.well_radius, ideal(params.depth), params.success_rate, params.cost_year)

    @staticmethod
    def waterBaseline(params = None, **kwargs):
        if params == None:
            params = SimulationParameters(**kwargs)
        return well(params.depth, params.well_radius, baseLine(params.depth), params.success_rate, params.cost_year)

    @staticmethod
    def waterIdeal(params = None, **kwargs):
        if params == None:
            params = SimulationParameters(**kwargs)
        return well(params.depth, params.well_radius, ideal(params.depth), params.success_rate, params.cost_year)

X_IC_well = 1.05
X_PC_well = 1.15

def baseLine(wellLength):
    well_typeA = 0.105*wellLength**2
    well_typeB = 1776.
    return (well_typeA, well_typeB)

def ideal(wellLength):
    well_typeA = -62.2*wellLength
    well_typeB = 1290.
    return (well_typeA, well_typeB)


def well(well_length, well_radius, well_type, success_rate, cost_year, dC_well = 0.):
    PPI_O_G = readCostTable('PPI_O&G', cost_year = cost_year)

    C_well = X_IC_well * X_PC_well * PPI_O_G * (well_type[0] + well_type[1] * 2 * well_radius * well_length + 275300.)
    return (C_well + dC_well) / success_rate

def wellCO2(well_length, well_radius, well_type, success_rate, cost_year):
    PPI_O_G = readCostTable('PPI_O&G', cost_year = cost_year)
    dC_well = X_IC_well * X_PC_well * PPI_O_G * (265. * 2 * well_radius * well_length + 133. * well_length)
    return well(well_length, well_radius, well_type, success_rate, cost_year, dC_well)
