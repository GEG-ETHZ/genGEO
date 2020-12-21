from enum import Enum

class WellFieldType(Enum):
    Doublet = 1
    _5Spot = 2
    _5Spot_SharedNeighbor = 3
    _5Spot_Many = 4

    def getInjWellMdotMultiplier(self):
        if self == WellFieldType.Doublet:
            return 1
        elif self == WellFieldType._5Spot:
            return 4
        else:
            return 4

    def getProdWellMdotMultiplier(self):
        if self == WellFieldType.Doublet:
            return 1
        elif self == WellFieldType._5Spot:
            return 1
        else:
            return 4