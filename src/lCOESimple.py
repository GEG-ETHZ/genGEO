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
import numpy as np

class LCOESimpleOutput(object):
    """LCOESimpleOutput."""
    pass

class LCOESimple(object):
    """LCOESimple."""
    def __init__(self, params = None, **kwargs):
        self.params = params
        if self.params == None:
            self.params = SimulationParameters(**kwargs)

        #Check heats & powers
        if self.params.F_OM < 0:
            raise Exception('GenGeo::LCOESimple:NegativeF_OM - Negative O&M Fraction!')
        elif self.params.discount_rate < 0:
            raise Exception('GenGeo::LCOESimple:NegativeDiscountRate - Negative Discount Rate!')
        elif self.params.lifetime <= 0:
            raise Exception('GenGeo::LCOESimple:Negativelifetime - Negative lifetime!')
        elif self.params.capacity_factor <= 0 or self.params.capacity_factor > 1:
            raise Exception('GenGeo::LCOESimple:Badcapacity_factor - Bad Capacity Factor!')

    def specificCapitalCost(self, capital_cost, capacity_W):
        if capacity_W <= 0 or capital_cost <= 0:
            return np.nan
        return capital_cost / capacity_W

    def lCOE(self, capital_cost, capacity_W):
        if capacity_W <= 0 or capital_cost <= 0:
            return np.nan
        CRF = self.params.discount_rate * (1 + self.params.discount_rate)**self.params.lifetime / ((1 + self.params.discount_rate)**self.params.lifetime - 1)
        return capital_cost * (CRF + self.params.F_OM) / (capacity_W * self.params.capacity_factor * 8760)

    def solve(self, capital_cost, energy_results):
        results = LCOESimpleOutput()
        capacity_W = energy_results.W_net_total
        results.LCOE = self.lCOE(capital_cost, capacity_W)
        results.specific_capital_cost = self.specificCapitalCost(capital_cost, capacity_W)
        self.output = results
        return results
