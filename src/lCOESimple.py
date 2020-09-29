import numpy as np

class LCOESimpleOutput(object):
    """LCOESimpleOutput."""
    pass

class LCOESimple(object):
    """LCOESimple."""
    def __init__(self, F_OM, discountRate, Lifetime, CapacityFactor):
        self.F_OM = F_OM
        self.discountRate = discountRate
        self.Lifetime = Lifetime
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

    def specificCapitalCost(self, CapitalCost):
        if self.Capacity_W <= 0 or CapitalCost <= 0:
            return np.nan
        return CapitalCost / self.Capacity_W

    def lCOE(self, CapitalCost):
        if self.Capacity_W <= 0 or CapitalCost <= 0:
            return np.nan
        CRF = self.discountRate * (1 + self.discountRate)**self.Lifetime / ((1 + self.discountRate)**self.Lifetime - 1)
        return CapitalCost * (CRF + self.F_OM) / (self.Capacity_W * self.CapacityFactor * 8760)

    def solve(self, CapitalCost, energy_results):
        self.Capacity_W = energy_results.W_net_total
        self.LCOE = self.lCOE(CapitalCost)
        self.specific_capital_cost = self.specificCapitalCost(CapitalCost)

    def gatherOutput(self):
        output = LCOESimpleOutput()
        output.LCOE = self.LCOE
        output.specific_capital_cost = self.specific_capital_cost
        return output
