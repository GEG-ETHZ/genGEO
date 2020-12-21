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

def testAssert(val1, val2, string, max_error = 1e-4):
    # relative error of the computed values compared to reference values provided by badams
    string += ' %.4e != %.4e'%(val1, val2)
    return (np.isclose([val1], [val2], rtol=max_error), string)
