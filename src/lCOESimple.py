import numpy as np

class LCOESimpleOutput(object):
    """LCOESimpleOutput."""
    pass

class LCOESimple(object):
    """LCOESimple."""
    def __init__(self, params = None, F_OM = None, discount_rate = None, lifetime = None, capacity_factor = None):
        if params:
            self.F_OM = params.F_OM
            self.discount_rate = params.discount_rate
            self.lifetime = params.lifetime
            self.capacity_factor = params.capacity_factor
        else:
            self.F_OM = F_OM
            self.discount_rate = discount_rate
            self.lifetime = lifetime
            self.capacity_factor = capacity_factor

        #Check heats & powers
        if self.F_OM < 0:
            raise ValueError('LCOE_Simple:NegativeF_OM - Negative O&M Fraction!')
        elif self.discount_rate < 0:
            raise ValueError('LCOE_Simple:NegativeDiscountRate - Negative Discount Rate!')
        elif self.lifetime <= 0:
            raise ValueError('LCOE_Simple:Negativelifetime - Negative lifetime!')
        elif self.capacity_factor <= 0 or self.capacity_factor > 1:
            raise ValueError('LCOE_Simple:Badcapacity_factor - Bad Capacity Factor!')

    def specificCapitalCost(self, CapitalCost):
        if self.Capacity_W <= 0 or CapitalCost <= 0:
            return np.nan
        return CapitalCost / self.Capacity_W

    def lCOE(self, CapitalCost):
        if self.Capacity_W <= 0 or CapitalCost <= 0:
            return np.nan
        CRF = self.discount_rate * (1 + self.discount_rate)**self.lifetime / ((1 + self.discount_rate)**self.lifetime - 1)
        return CapitalCost * (CRF + self.F_OM) / (self.Capacity_W * self.capacity_factor * 8760)

    def solve(self, CapitalCost, energy_results):
        self.Capacity_W = energy_results.W_net_total
        self.LCOE = self.lCOE(CapitalCost)
        self.specific_capital_cost = self.specificCapitalCost(CapitalCost)

    def gatherOutput(self):
        output = LCOESimpleOutput()
        output.LCOE = self.LCOE
        output.specific_capital_cost = self.specific_capital_cost
        return output
