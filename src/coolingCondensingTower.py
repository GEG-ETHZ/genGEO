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

from utils.constantsAndPaths import ConversionConstants
import numpy as np

from models.coolingCondensingTowerMode import CoolingCondensingTowerMode

from utils.readXlsxData import readCostTable


class CoolingCondensingTower(object):   
    """CoolingCondensingTower."""

    @staticmethod
    def computeProcess(a, b, c, d, T_ambient_C, dT_approach_CT, dT_range_CT):
        A = a * (1/dT_approach_CT)
        B = b * (T_ambient_C+ConversionConstants.kelvin2celsius)
        C = c * (T_ambient_C+ConversionConstants.kelvin2celsius)/dT_approach_CT
        D = d * (1/(dT_approach_CT + dT_range_CT))
        return A + B + C + D

    @staticmethod
    def specificCaptitalCost(Q_cooler, Q_condenser, TDC, T_ambient_C, dT_approach_CT, dT_range_CT, cost_year, coolingMode):
        def processWet(process):
            if process == 'cooling':
                return CoolingCondensingTower.computeProcess(5.58e3, 0., -1.77e1, 1.96e2, T_ambient_C, dT_approach_CT, dT_range_CT)
            elif process == 'condensing':
                return CoolingCondensingTower.computeProcess(4.08e3, -1.54e-2, -1.24e1, 0., T_ambient_C, dT_approach_CT, dT_range_CT)
            else:
                raise Exception('GenGeo::coolingCondensingTower:UnknownProcess - Unknown Process - use "cooling" or "condensing"')

        def processDry(process):
            if process == 'cooling':
                return CoolingCondensingTower.computeProcess(7.31e3, 0., 0., 1.23e3, T_ambient_C, dT_approach_CT, dT_range_CT)
            elif process == 'condensing':
                return CoolingCondensingTower.computeProcess(1.91e3, 0., 0., 0., T_ambient_C, dT_approach_CT, dT_range_CT)
            else:
                raise Exception('GenGeo::coolingCondensingTower:UnknownProcess - Unknown Process - use "cooling" or "condensing"')

        if coolingMode == CoolingCondensingTowerMode.Wet:
            c_cooling = processWet('cooling')
            c_condensing = processWet('condensing')
        elif coolingMode == CoolingCondensingTowerMode.Dry:
            c_cooling = processDry('cooling')
            c_condensing = processDry('condensing')
        else:
            raise Exception('GenGeo::coolingCondensingTower:UnknownCoolingMode - Unknown Cooling Mode')

        # c_cooling and c_condensing are both in units of $/kWth
        if np.isnan(Q_cooler):
            Q_cooler = 0
        if np.isnan(Q_condenser):
            Q_condenser = 0

        # Reference case 1000 kWth (1e6 Wth)
        Q_Ref_BAC = 1e6
        F_cooling = abs(Q_cooler)/(abs(Q_cooler)+abs(Q_condenser))
        C_Ref_BAC = readCostTable('PPI_PE', cost_year = cost_year) * Q_Ref_BAC * TDC * (F_cooling*(c_cooling/1e3) + (1-F_cooling)*(c_condensing/1e3))
        return C_Ref_BAC * (abs(Q_cooler+Q_condenser)/Q_Ref_BAC)**0.8

    @staticmethod
    def parasiticPowerFraction(T_ambient_C, dT_approach_CT, dT_range_CT, coolingMode):
        def processWet(process):
            if process == 'cooling':
                return CoolingCondensingTower.computeProcess(1.2, 0., -3.79e-3, 1.95e-2, T_ambient_C, dT_approach_CT, dT_range_CT)
            elif process == 'condensing':
                return CoolingCondensingTower.computeProcess(1.65, -6.24e-6, -5.03e-3, 0., T_ambient_C, dT_approach_CT, dT_range_CT)
            else:
                raise Exception('GenGeo::coolingCondensingTower:UnknownProcess - Unknown Process - use "cooling" or "condensing"')

        def processDry(process):
            if process == 'cooling':
                return CoolingCondensingTower.computeProcess(7.65e-1, 0., 0., 1.28e-1, T_ambient_C, dT_approach_CT, dT_range_CT)
            elif process == 'condensing':
                return CoolingCondensingTower.computeProcess(0.619, 0., 0., 0., T_ambient_C, dT_approach_CT, dT_range_CT)
            else:
                raise Exception('GenGeo::coolingCondensingTower:UnknownProcess - Unknown Process - use "cooling" or "condensing"')

        if coolingMode == CoolingCondensingTowerMode.Wet:
            return processWet
        elif coolingMode == CoolingCondensingTowerMode.Dry:
            return processDry
        else:
            raise Exception('GenGeo::coolingCondensingTower:UnknownCoolingMode - Unknown Cooling Mode')