from pathlib import Path
import os

def getProjectRoot() -> Path:
    return Path(__file__).parent.parent

def getWellCost():
    return os.path.join(getProjectRoot(), 'data', 'PPI_Table.xlsx')

def getTboilOptimum(orcFluid):
    return os.path.join(getProjectRoot(), 'data', 'ORC_Tboil_optimum_%s.csv'%orcFluid)

class globalConstants(object):
    """globalConstants carries global constants to convert units."""

    secPerYear = 3600. * 24. * 365.25
    kelvin2celsius = 273.15
