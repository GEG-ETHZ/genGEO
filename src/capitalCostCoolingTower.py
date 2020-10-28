
import numpy as np

from utils.readXlsxData import readCostTable


class CapitalCostCoolingTower(object):
    """CapitalCostCoolingTower."""
    @staticmethod
    def cost(Q_cooler, Q_condenser, TDC, T_ambient_C, dT_approach_CT, dT_range_CT, cost_year):
        a_cool = 5.58e3
        b_cool = 0
        c_cool = -1.77e1
        d_cool = 1.96e2
        c_cooling_wet = a_cool*(1/dT_approach_CT) + b_cool*(T_ambient_C+273.15) + c_cool*(T_ambient_C+273.15)/dT_approach_CT + d_cool*(1/(dT_approach_CT + dT_range_CT))
        a_cond = 4.08e3
        b_cond = -1.54e-2
        c_cond = -1.24e1
        d_cond = 0
        c_condensing_wet = a_cond*(1/dT_approach_CT) + b_cond*(T_ambient_C+273.15) + c_cond*(T_ambient_C+273.15)/dT_approach_CT + d_cond*(1/(dT_approach_CT + dT_range_CT))

        # c_cooling_wet and c_condensing_wet both in units of $/kWth

        # Reference case 1000 kWth (1e6 Wth)
        Q_Ref_BAC = 1e6
        F_cooling = abs(Q_cooler)/(abs(Q_cooler)+abs(Q_condenser))
        C_Ref_BAC = readCostTable('PPI_PE', cost_year = cost_year) * Q_Ref_BAC * TDC * (F_cooling*(c_cooling_wet/1e3) + (1-F_cooling)*(c_condensing_wet/1e3))
        return C_Ref_BAC * (abs(Q_cooler+Q_condenser)/Q_Ref_BAC)**0.8
