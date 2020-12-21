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

from utils.fluidState import FluidState

class PorousReservoirBase(object):
    """PorousReservoirBase."""

    def __init__(self):
        self.params.T_reservoir = lambda : abs(self.params.depth) * abs(self.params.dT_dz) + self.params.T_surface_rock
        self.params.P_reservoir = lambda : abs(self.params.depth) * 1000. * self.params.g
        self.params.P_reservoir_max = lambda : abs(self.params.depth) * 2500. * self.params.g
        self.params.P_system_min = lambda : FluidState.getStateFromTQ(self.params.T_reservoir(), 0, self.params.working_fluid).P_Pa + 1e5
