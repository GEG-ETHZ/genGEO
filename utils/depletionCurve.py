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
import math

def depletionCurve(psi, p1, p2, p3):
    if psi <= (-p1/p2):
        return 1 - (1 - p3) * (1 - 0.5*math.erfc(p2*psi + p1))
    else:
        C = 1.13
        return 1 - (1 - p3) * (1 - 0.5*np.exp(-C*(p2*psi + p1)))
