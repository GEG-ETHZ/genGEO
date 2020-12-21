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


class PowerPlantEnergyOutput(object):
    """PowerPlantEnergyOutput."""
    q_preheater = np.nan
    q_recuperator = np.nan
    q_boiler = np.nan
    q_desuperheater = np.nan
    q_condenser = np.nan
    w_turbine = np.nan
    w_pump = np.nan
    w_cooler = np.nan
    w_condenser = np.nan
    w_net = np.nan


class PowerPlantOutput(PowerPlantEnergyOutput):
    """PowerPlantOutput."""
    dT_range_CT = np.nan
    dT_LMTD_preheater = np.nan
    dT_LMTD_recuperator = np.nan
    dT_LMTD_boiler = np.nan
    state = np.nan
