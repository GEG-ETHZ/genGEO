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
from pathlib import Path
import os
import numpy as np

from models.optimizationType import OptimizationType

def getProjectRoot() -> Path:
    return Path(__file__).parent.parent

def getWellCost():
    return os.path.join(getProjectRoot(), 'data', 'PPI_Table.xlsx')

def getPboilOptimum():
    return os.path.join(getProjectRoot(), 'data', 'ORC_Pboil_optimum.csv')

def getTboilOptimum():
    path = os.path.join(getProjectRoot(), 'data', 'ORC_Tboil_optimum_%s_%s.csv')
    data_dict = {}
    for opt_mod in [ OptimizationType.MaximizePower, OptimizationType.MinimizeCost ]:
        data_dict[opt_mod] = {}
        for orc_fluid in ['R600a', 'R245fa']:
            if opt_mod == OptimizationType.MaximizePower:
                data_dict[opt_mod][orc_fluid] = np.genfromtxt(path%('maxPower', orc_fluid), delimiter=',')
            elif opt_mod == OptimizationType.MinimizeCost:
                data_dict[opt_mod][orc_fluid] = np.genfromtxt(path%('minCost', orc_fluid), delimiter=',')
            
    return data_dict

class ConversionConstants(object):
    """ConversionConstants carries global constants."""

    secPerYear = 3600. * 24. * 365.
    kelvin2celsius = 273.15
