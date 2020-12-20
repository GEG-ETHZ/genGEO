import numpy as np

from utils.fluidState import FluidState

def maxSubcritORCBoilTemp(orc_fluid):
    P_crit = FluidState.getPcrit(orc_fluid)

    # assume min pressure is condensing at 0C
    P_min = FluidState.getStateFromTQ(0, 1, orc_fluid).P_Pa
    dP = (P_crit-P_min)/1000

    # find pressure with maximum entropy
    P = np.arange(P_min, P_crit, dP)
    s = FluidState.getStateFromPQ(P, 1, orc_fluid).s_JK

    max_s = np.max(s)
    row_index = np.argmax(s, axis=0)
    max_P = P[row_index]

    return FluidState.getStateFromPS(max_P, max_s, orc_fluid).T_C
