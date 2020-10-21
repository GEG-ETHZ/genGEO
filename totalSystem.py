import numpy as np

from src.totalSystemWater import TotalSystemWater
from src.totalSystemCO2 import TotalSystemCO2

import time

logTrans = np.arange(2., 8., 1.)
permeabilities = 1e-15 * 10. ** logTrans
depths = np.arange(1000, 8000, 1000)

permeabilities = [1e-9]
depths = [5000]
file = open('data_tmp.csv', 'w')

# file = open('data_Water_R600a_orig.csv', 'w')
# file = open('data_Water_R245fa_orig.csv', 'w')
# file = open('data_CO2_orig.csv', 'w')

for depth in depths:
    for permeability in permeabilities:
        start = time.time()
        print('######### new case ########')
        print(depth)
        print(permeability)
        # system = TotalSystemWater(depth = depth, permeability = permeability / 100., orc_fluid = 'R600a')
        system = TotalSystemWater(depth = depth, permeability = permeability / 100., orc_fluid = 'R245fa')
        # system = TotalSystemCO2(depth = depth, permeability = permeability / 100.)
        try:
            results = system.minimizeLCOEBrownfield(time_years = 1.)
            error_str =  'No error'
            lcoe_b = results.capital_cost_model.LCOE_brownfield.LCOE * 1e6
            lcoe_g = results.capital_cost_model.LCOE_greenfield.LCOE * 1e6
            power = results.energy_results.W_net / 1e6
            optMdot = results.optMdot

        except Exception as error:
            error_str = str(error).replace("\n", "").replace(",", " - ")
            lcoe_b = np.nan
            lcoe_g = np.nan
            power = np.nan
            optMdot = np.nan

        print(error_str)
        print('optMdot ', optMdot)
        print('lcoe_b ', lcoe_b)
        print('lcoe_g ', lcoe_g)
        print('power ', power)

        if np.isnan(lcoe_b):
            lcoe_b = 0.
        if np.isnan(lcoe_g):
            lcoe_g = 0.
        if np.isnan(optMdot):
            optMdot = 0.
        if np.isnan(power):
            power = 0.


        end = time.time()
        print('time', end - start)
        file.write(','.join([str(i) for i in [depth, permeability, optMdot, lcoe_b, lcoe_g, power, """%s\n"""%error_str]]))

file.close()
