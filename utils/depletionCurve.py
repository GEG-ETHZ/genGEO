import numpy as np
import math

def depletionCurve(psi, p1, p2, p3):
    if psi <= (-p1/p2):
        return 1 - (1 - p3) * (1 - 0.5*math.erfc(p2*psi + p1))
    else:
        C = 1.13
        return 1 - (1 - p3) * (1 - 0.5*np.exp(-C*(p2*psi + p1)))
