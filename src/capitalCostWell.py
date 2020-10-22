from utils.readXlsxData import readCostTable


class CapitalCostWell(object):
    """CapitalCostWell."""
    @staticmethod
    def cO2Baseline(well_length = None, well_radius = None, success_rate = None, cost_year = None, params = None):
        if params:
            well_length = params.well_length
            well_radius = params.well_radius
            success_rate = params.success_rate
            cost_year = params.cost_year
        return wellCO2(well_length, well_radius, baseLine(well_length), success_rate, cost_year)

    @staticmethod
    def cO2Ideal(well_length = None, well_radius = None, success_rate = None, cost_year = None, params = None):
        if params:
            well_length = params.well_length
            well_radius = params.well_radius
            success_rate = params.success_rate
            cost_year = params.cost_year
        return wellCO2(well_length, well_radius, ideal(well_length), success_rate, cost_year)

    @staticmethod
    def waterBaseline(well_length = None, well_radius = None, success_rate = None, cost_year = None, params = None):
        if params:
            well_length = params.well_length
            well_radius = params.well_radius
            success_rate = params.success_rate
            cost_year = params.cost_year
        return well(well_length, well_radius, baseLine(well_length), success_rate, cost_year)

    @staticmethod
    def waterIdeal(well_length = None, well_radius = None, success_rate = None, cost_year = None, params = None):
        if params:
            well_length = params.well_length
            well_radius = params.well_radius
            success_rate = params.success_rate
            cost_year = params.cost_year
        return well(well_length, well_radius, ideal(well_length), success_rate, cost_year)

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
    PPI_O_G = readCostTable(cost_year, 'PPI_O&G')

    C_well = X_IC_well * X_PC_well * PPI_O_G * (well_type[0] + well_type[1] * 2 * well_radius * well_length + 275300.)
    return (C_well + dC_well) / success_rate

def wellCO2(well_length, well_radius, well_type, success_rate, cost_year):
    PPI_O_G = readCostTable(cost_year, 'PPI_O&G')
    dC_well = X_IC_well * X_PC_well * PPI_O_G * (265. * 2 * well_radius * well_length + 133. * well_length)
    return well(well_length, well_radius, well_type, success_rate, cost_year, dC_well)
