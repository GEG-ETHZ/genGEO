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
from enum import Enum

class WellFieldType(Enum):
    Doublet = 1
    _5Spot = 2
    _5Spot_SharedNeighbor = 3
    _5Spot_Many = 4

    def getInjWellMdotMultiplier(self):
        if self == WellFieldType.Doublet:
            return 1
        elif self == WellFieldType._5Spot:
            return 4
        else:
            return 4

    def getProdWellMdotMultiplier(self):
        if self == WellFieldType.Doublet:
            return 1
        elif self == WellFieldType._5Spot:
            return 1
        else:
            return 4