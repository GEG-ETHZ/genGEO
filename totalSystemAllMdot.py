import numpy as np

from src.totalSystemWater import TotalSystemWater
from src.totalSystemCO2 import TotalSystemCO2

import time

logTrans = np.arange(2., 8., 1.)
permeabilities = 1e-15 * 10. ** logTrans
depths = np.arange(1000, 8000, 1000)

for depth in depths:
    for permeability in permeabilities:
        start = time.time()
        filename = 'data_tmp_test_%s_%.e.csv'%(depth, permeability)
        print(filename)
        file = open(filename, 'w')
        for m_dot in range(1, 200):
            print('\n### new m_dot %s ###'%m_dot)
            system = TotalSystemWater(depth = depth, permeability = permeability / 100., orc_fluid = 'R245fa')
            # system = TotalSystemWater(depth = 1000, permeability = 1e-13 / 100.)
            # system = TotalSystemCO2(depth = 1000, permeability = 1e-12 / 100.)
            try:
                results = system.solve(m_dot = m_dot, time_years = 1.)
                lcoe = results.capital_cost_model.LCOE_brownfield.LCOE * 1e6

                lcoe_b = results.capital_cost_model.LCOE_brownfield.LCOE * 1e6
                lcoe_g = results.capital_cost_model.LCOE_greenfield.LCOE * 1e6
                power = results.energy_results.W_net / 1e6
            except Exception as ex:
                error_str = str(ex)
                lcoe_b = np.nan
                lcoe_g = np.nan
                power = np.nan

            file.write(','.join([str(i) for i in [m_dot, lcoe_b, lcoe_g, power, """%s\n"""%error_str]]))

            print(error_str)
            print('lcoe_b ', lcoe_b)
            print('lcoe_g ', lcoe_g)
            print('power ', power)

        file.close()

        end = time.time()
        print('time', end - start)
