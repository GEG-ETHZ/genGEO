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

from utils.getIntersection import getIntersection

class Solver(object):
    """Solver."""

    def __init__(self):
        self.func1 = []
        self.allowExtrapolation = True
        self.showPlot = False

    def addDataAndEstimate(self, x, y):
        self.func1.append(np.array([x, y]))
        return self.convergeToZero()

    def convergeToZero(self):

        func1 = np.array(self.func1)
        func1 = func1[np.argsort(func1[:, 0])]
        func1 = np.unique(func1, axis=0)
        func2 = np.array([np.array([0, 0]), np.array([1, 0])])

        try:
            x_zero, y_zero = getIntersection(func1, func2, self.allowExtrapolation)
        except Exception as ex:
            if str(ex).find('GetIntersection:NoIntersection') > -1:
                # The lines appear divergent on either end.
                x_zero = np.nan
                y_zero = 0
                print('Converge to zero not working. Nudging x to %s' %x_zero)
            elif str(ex).find('GetIntersection:NotEnoughRows') > -1:
                # Only single values
                # Guess last result
                x_zero = np.nan
                if func1[-1, 1] == 0:
                    x_zero = func1[-1, 0]
                y_zero = 0
            else:
                raise ex
        return x_zero
