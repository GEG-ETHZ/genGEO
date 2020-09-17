import numpy as np


from utils.readXlsxData import readCostTable

class CapitalCostSurfacePlantCPGResults(object):
    """CapitalCostSurfacePlantCPGResults."""
    def __init__(self):


class CapitalCostSurfacePlantCPG(object):
    """CapitalCostSurfacePlantCPG."""

    def __init__(self, cost_year):
        self.cost_year = cost_year
        self.ppi_T_G = readCostTable(cost_year, 'PPI_T-G')
        self.ppi_pump = readCostTable(cost_year, 'PPI_Pump')

    def solve(W_turbine, Q_desuperheater, Q_condenser, W_pump_inj, T_ambient_C,
                dT_approach_CT, dT_range_CT):

        # # TODO: Do some checks if temp is correct

        results = CapitalCostSurfacePlantCPGResults()

        if 0. < W_pump_inj < 10000.: W_pump_inj = 0.

        # C_T_G
        # Regular fluid
        S_T_fluid = 1.20 #CO2
        results.C_T_G = 0.67 * PPI_T_G * (S_T_fluid*2830*(W_turbine/1e3)**0.745 + 3680*(W_turbine/1e3)**0.617)

        # C_pump_inj
        C_pump_surface_inj = 1750 * (1.34*-1*(W_pump_inj/1e3))**0.7
        S_pump_inj = 2.09 #CO2
        results.C_pump_inj = PPI_pump * S_pump_inj * C_pump_surface_inj

        # C_coolingTowers
        TDC = 1.2
        results.C_coolingTowers= CapitalCostCoolingTower.cost(Q_desuperheater, Q_condenser, TDC,
                                                    T_ambient_C, dT_approach_CT, dT_range_CT, self.cost_year)

        # C_primaryEquipment
        C_primaryEquipment = results.C_T_G + results.C_pump_inj + results.C_coolingTowers

        X_SE = 1.39
        X_CL = 0.58
        X_CM = 0.11
        X_ST = 0.00
        X_F = 0.04
        X_PC_sp = 1.15
        X_IC_sp = 1.12

        C_plant_TEC = X_SE * C_primaryEquipment
        results.C_plant_otherEquipment = C_plant_TEC - C_primaryEquipment

        C_plant_BEC = C_plant_TEC * (1 + X_CL + X_CM + X_ST + X_F)
        results.C_plant_installation = C_plant_BEC - C_plant_TEC

        results.C_plant = X_PC_sp * X_IC_sp * C_plant_BEC
        results.C_plant_indirectContingency = results.C_plant - C_plant_BEC

        # calculate fractions
        results.f_T_G = results.C_T_G / C_primaryEquipment
        results.f_pump_inj = results.C_pump_inj / C_primaryEquipment
        results.f_coolingTowers = results.C_coolingTowers / C_primaryEquipment
