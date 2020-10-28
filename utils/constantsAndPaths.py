from pathlib import Path
import os
import numpy as np


def getProjectRoot() -> Path:
    return Path(__file__).parent.parent

def getWellCost():
    return os.path.join(getProjectRoot(), 'data', 'PPI_Table.xlsx')

def getPboilOptimum():
    return os.path.join(getProjectRoot(), 'data', 'ORC_Pboil_optimum.csv')

def getTboilOptimum():
    path = os.path.join(getProjectRoot(), 'data', 'ORC_Tboil_optimum_%s.csv')
    R600a = np.genfromtxt(path%'R600a', delimiter=',')
    R245fa = np.genfromtxt(path%'R245fa', delimiter=',')
    return {'R600a': R600a, 'R245fa':R245fa}

class ConversionConstants(object):
    """ConversionConstants carries global constants."""

    secPerYear = 3600. * 24. * 365.
    kelvin2celsius = 273.15
