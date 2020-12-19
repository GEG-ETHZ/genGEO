
from utils.readXlsxData import readCostTable
from utils.simulationParameters import SimulationParameters

class CapitalCostSurfacePipes(object):
    """CapitalCostSurfacePipes."""

    @staticmethod
    def cost(params = None, **kwargs):
        if params == None:
            params = SimulationParameters(**kwargs)

        X_PCs = 1.15
        X_ICs = 1.12

        if params.wellFieldType == 'Doublet':
            L_surfacePipe = 707
            D_surfacePipe = 0.41
        elif params.wellFieldType == '5spot':
            L_surfacePipe = 3000
            D_surfacePipe = 0.41
        elif params.wellFieldType == '5spot_SharedNeighbor':
            L_surfacePipe = 707
            D_surfacePipe = 0.41
        elif params.wellFieldType == '5spot_ManyN':
            L_surfacePipe_manyN = {
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
            D_surfacePipe_manyN = {
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
            L_surfacePipe = L_surfacePipe_manyN[params.N_5spot]
            D_surfacePipe = D_surfacePipe_manyN[params.N_5spot]
        else:
            L_surfacePipe = 707
            D_surfacePipe = 0.41

        c_surfacePipe = 2205. * D_surfacePipe**2 + 134.
        return X_PCs * X_ICs * readCostTable('PPI_Pipe', cost_year = params.cost_year) * c_surfacePipe * L_surfacePipe
