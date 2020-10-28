import pandas as pd

from utils.constantsAndPaths import getWellCost

def readXlsxColumn(file, sheet, headers, headerline):
    return pd.read_excel(file, sheet_name = sheet, header = headerline, usecols = headers)

def readCostTable(column, cost_year = None):
    table = readXlsxColumn(getWellCost(), 'Sheet1', ['Year', column], 1)
    if cost_year == None:
        return dict(table.to_numpy())
    return table.query('Year == %s'%cost_year)[column].values[0]
