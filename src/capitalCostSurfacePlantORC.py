
from utils.readXlsxData import readCostTable

class CapitalCostSurfacePlantORCResults(object):
    """CapitalCostSurfacePlantORCResults."""
    def __init__(self):


class CapitalCostSurfacePlantORC(object):
    """CapitalCostSurfacePlantORC."""

    def __init__(self, cost_year):
        self.cost_year = cost_year
        self.ppi_T_G = readCostTable(cost_year, 'PPI_T-G')
        self.ppi_pump = readCostTable(cost_year, 'PPI_Pump')
        self.ppi_HX = readCostTable(cost_year, 'PPI_HX')

    def solve(Q_preheater, Q_boiler, W_turbine, Q_recuperator, Q_desuperheater, Q_condenser,
                W_pump_orc, W_pump_prod, T_ambient_C,
                dT_approach_CT, dT_range_CT, dT_LMTD_preheater, dT_LMTD_boiler, dT_LMTD_recuperator):
                
        # # TODO: Do some checks if temp is correct

        results = CapitalCostSurfacePlantORCResults()

        # C_T_G
        #Regular fluid
        S_T_fluid = 1.00  #Not CO2
        results.C_T_G = 0.67 * self.ppi_T_G * (S_T_fluid*2830*(W_turbine/1e3)**0.745 + 3680*(W_turbine/1e3)**0.617)

        # C_pump (ORC)
        C_pump_orc_surface = 1750 * (1.34*-1*(W_pump_orc/1e3))**0.7
        S_pump_orc = 1.00  #Water
        results.C_pump_orc = self.ppi_pump * S_pump_orc * C_pump_orc_surface

        # C_coolingTowers
        TDC = 1
        results.C_coolingTowers = CapitalCostCoolingTower.cost(Q_desuperheater, Q_condenser, TDC,
                                                    T_ambient_C, dT_approach_CT, dT_range_CT, self.cost_year)

        # C_heatExchanger
        #dT_LMTD_HX
        #U = 500/1000 #kW/m**2-K
        U = 500 #W/m**2-K
        if (isnan(dT_LMTD_preheater) || dT_LMTD_preheater == 0)
            A_preheater = 0
        else
            A_preheater = Q_preheater / U / dT_LMTD_preheater
        end
        A_boiler = Q_boiler / U / dT_LMTD_boiler
        A_HX = A_preheater + A_boiler
        results.C_heatExchanger = self.ppi_HX * (239*A_HX + 13400)

        # C_recuperator
        if (isnan(dT_LMTD_recuperator) || dT_LMTD_recuperator == 0)
            A_recuperator = 0
        else
            A_recuperator = Q_recuperator / U / dT_LMTD_recuperator
        end
        results.C_recuperator = self.ppi_HX * (239*A_recuperator + 13400)

        # C_productionPump
        C_pump_prod_lineshaft = 1750 * (1.34*-1*(W_pump_prod/1e3))**0.7 + 5750 * (1.34*-1*(W_pump_prod/1e3))**0.2
        S_pump_prod = 1.00 #Water
        results.C_pump_prod = self.ppi_pump * S_pump_prod * C_pump_prod_lineshaft

        # THEN
        # C_primaryEquipment
        C_primaryEquipment = results.C_T_G + results.C_pump_orc + results.C_coolingTowers + results.C_heatExchanger + results.C_pump_prod + results.C_recuperator

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
        results.C_plant_otherEquipment = C_plant_BEC - C_plant_TEC

        results.C_plant = X_PC_sp * X_IC_sp * C_plant_BEC
        results.C_plant_indirectContingency = results.C_plant - C_plant_BEC

        # calculate fractions
        results.f_T_G = results.C_T_G / C_primaryEquipment
        results.f_pump_orc = results.C_pump_orc / C_primaryEquipment
        results.f_coolingTowers = results.C_coolingTowers / C_primaryEquipment
        results.f_heatExchanger = results.C_heatExchanger / C_primaryEquipment
        results.f_pump_prod = results.C_pump_prod / C_primaryEquipment

        return results
