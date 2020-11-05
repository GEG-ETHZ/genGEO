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
    path = os.path.join(getProjectRoot(), 'data', 'ORC_Tboil_optimum_%s_%s.csv')
    data_dict = {}
    for opt_mod in ['maxPower', 'minCost']:
        data_dict[opt_mod] = {}
        for orc_fluid in ['R600a', 'R245fa']:
            data_dict[opt_mod][orc_fluid] = np.genfromtxt(path%(opt_mod, orc_fluid), delimiter=',')
    return data_dict

class ConversionConstants(object):
    """ConversionConstants carries global constants."""

    secPerYear = 3600. * 24. * 365.
    kelvin2celsius = 273.15
