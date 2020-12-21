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
import pandas as pd

from utils.constantsAndPaths import getWellCost

def readXlsxColumn(file, sheet, headers, headerline):
    return pd.read_excel(file, sheet_name = sheet, header = headerline, usecols = headers)

def readCostTable(column, cost_year = None):
    table = readXlsxColumn(getWellCost(), 'Sheet1', ['Year', column], 1)
    if cost_year == None:
        return dict(table.to_numpy())
    return table.query('Year == %s'%cost_year)[column].values[0]
