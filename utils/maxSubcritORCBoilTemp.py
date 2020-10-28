import numpy as np

from utils.fluidStates import FluidState

def maxSubcritORCBoilTemp(orc_fluid):
    P_crit = FluidState.getPcrit(orc_fluid)

    # assume min pressure is condensing at 0C
    P_min = FluidState.getPFromTQ(0, 1, orc_fluid)
    dP = (P_crit-P_min)/1000

    # find pressure with maximum entropy
    P = np.arange(P_min, P_crit, dP)
    s = FluidState.getSFromPQ(P, 1, orc_fluid)

    max_s = np.max(s)
    row_index = np.argmax(s, axis=0)
    max_P = P[row_index]

    return FluidState.getTFromPS(max_P, max_s, orc_fluid)
