import numpy as np

class LCOESimple(object):
    """LCOESimple."""
    def __init__(self, CapitalCost, F_OM, discountRate, Lifetime, Capacity_W, CapacityFactor):
        self.CapitalCost = CapitalCost
        self.F_OM = F_OM
        self.discountRate = discountRate
        self.Lifetime = Lifetime
        self.Capacity_W = Capacity_W
        self.CapacityFactor = CapacityFactor

        #Check heats & powers
        if F_OM < 0:
            raise ValueError('LCOE_Simple:NegativeF_OM - Negative O&M Fraction!')
        elif discountRate < 0:
            raise ValueError('LCOE_Simple:NegativeDiscountRate - Negative Discount Rate!')
        elif Lifetime <= 0:
            raise ValueError('LCOE_Simple:NegativeLifetime - Negative Lifetime!')
        elif CapacityFactor <= 0 or CapacityFactor > 1:
            raise ValueError('LCOE_Simple:BadCapacityFactor - Bad Capacity Factor!')

    def specificCapitalCost(self):
        if self.Capacity_W <= 0 or self.CapitalCost <= 0:
            return np.nan
        return self.CapitalCost / self.Capacity_W

    def lCOE(self):
        if self.Capacity_W <= 0 or self.CapitalCost <= 0:
            return np.nan
        CRF = self.discountRate * (1 + self.discountRate)**self.Lifetime / ((1 + self.discountRate)**self.Lifetime - 1)
        return self.CapitalCost * (CRF + self.F_OM) / (self.Capacity_W * self.CapacityFactor * 8760)
