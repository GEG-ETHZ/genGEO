
from utils.readXlsxData import readCostTable

class CapitalCostSurfacePipes(object):
    """CapitalCostSurfacePipes."""

    def __init__(self, N):
        self.N = N

    def solve(self, cost_year):
        X_PCs = 1.15
        X_ICs = 1.12

        L_surfacePipe = {-1: 707,
                        0: 3000,
                        1: 3000,
                        2: 12000,
                        3: 25000,
                        4: 45000,
                        5: 69000,
                        6: 107000,
                        7: 153000,
                        8: 223000,
                        9: 309000,
                        10: 406000
                        }

        D_surfacePipe = {-1: 0.41,
                        0: 0.34,
                        1: 0.34,
                        2: 0.54,
                        3: 0.65,
                        4: 0.79,
                        5: 0.89,
                        6: 0.98,
                        7: 1.02,
                        8: 1.05,
                        9: 1.09,
                        10: 1.13
                        }
        c_surfacePipe = 2205. * D_surfacePipe[self.N]**2 + 134.
        return X_PCs * X_ICs * readCostTable(cost_year, 'PPI_Pipe') * c_surfacePipe * L_surfacePipe[self.N]
