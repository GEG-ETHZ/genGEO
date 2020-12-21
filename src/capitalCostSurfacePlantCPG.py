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

from src.coolingCondensingTower import CoolingCondensingTower

from models.simulationParameters import SimulationParameters
from utils.readXlsxData import readCostTable

class CapitalCostSurfacePlantCPGResults(object):
    """CapitalCostSurfacePlantCPGResults."""
    pass


class CapitalCostSurfacePlantCPG(object):
    """CapitalCostSurfacePlantCPG."""

    def __init__(self, params = None, **kwargs):
        self.params = params
        if self.params == None:
            self.params = SimulationParameters(**kwargs)
        self.ppi_T_G = readCostTable('PPI_T-G')
        self.ppi_pump = readCostTable('PPI_Pump')

    def solve(self, energy_results, fluid_system):

        dT_range_CT =  0.
        T_ambient_C = self.params.T_ambient_C
        dT_approach_CT = self.params.dT_approach

        Q_condenser = energy_results.Q_condenser_total
        Q_desuperheater = energy_results.Q_desuperheater_total
        W_pump_inj = energy_results.W_pump_total
        W_turbine = energy_results.W_turbine_total

        if 0. < W_pump_inj < 10000.: W_pump_inj = 0.

        # Check heats & powers
        if W_turbine < 0:
            raise Exception('GenGeo::CapitalCostSurfacePlantCPG:NegativeTurbinePower - Negative Turbine Power')
        elif Q_desuperheater > 0 or Q_condenser > 0:
            raise Exception('GenGeo::CapitalCostSurfacePlantCPG:PositiveCondenserHeat - Positive Condenser Heat')
        elif W_pump_inj > 0:
            raise Exception('GenGeo::CapitalCostSurfacePlantCPG:PositiveCPGPumpPower - Positive CPG Pump Power')

        # Check temps
        if dT_approach_CT < 0 or dT_range_CT < 0:
            raise Exception('GenGeo::CapitalCostSurfacePlantCPG:NegativedT - Negative Temp Difference')

        results = CapitalCostSurfacePlantCPGResults()

        # C_T_G
        # Regular fluid
        S_T_fluid = 1.20 #CO2
        results.C_T_G = 0.67 * self.ppi_T_G[self.params.cost_year] * (S_T_fluid*2830*(W_turbine/1e3)**0.745 + 3680*(W_turbine/1e3)**0.617)

        # C_pump_inj
        C_pump_surface_inj = 1750 * (1.34*-1*(W_pump_inj/1e3))**0.7
        S_pump_inj = 2.09 #CO2
        results.C_pump_inj = self.ppi_pump[self.params.cost_year] * S_pump_inj * C_pump_surface_inj

        # C_coolingTowers
        TDC = 1.2
        results.C_coolingTowers = CoolingCondensingTower.specificCaptitalCost(Q_desuperheater, Q_condenser, TDC,
                                                    T_ambient_C, dT_approach_CT, dT_range_CT, self.params.cost_year, self.params.cooling_mode)

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

        return results
