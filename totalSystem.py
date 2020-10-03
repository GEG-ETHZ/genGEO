import numpy as np

from src.totalSystemWater import TotalSystemWater

class Dummy(object):
    pass


# m_dot = 0.768538599999999
# system = TotalSystemWater(depth = 5000, permeability = 1e-13 / 100.)
system = TotalSystemWater(depth = 7000, permeability = 1e-11 / 100.)
results = system.minimizeLCOEBrownfield(time_years = 1.)
# results = system.solve(m_dot = m_dot, time_years = 1.)
# results.optMdot = m_dot
print(results.optMdot)
print(results.capital_cost_model.LCOE_brownfield.LCOE * 1e6)
1/0

logTrans = np.arange(2., 8., 1.)
permeabilities = 1e-15 * 10. ** logTrans
depths = np.arange(1000, 8000, 1000)

file = open('test_out.csv', 'w')

for depth in depths:
    for permeability in permeabilities:
        print('######### new case ########')
        print(depth)
        print(permeability)
        system = TotalSystemWater(depth = depth, permeability = permeability / 100.)
        try:
            results = system.minimizeLCOEBrownfield(time_years = 1.)
            error_str =  'No error'
        except ValueError as error:
            error_str = str(error).replace("\n", "").replace(",", " - ")
            results = Dummy()
            results.optMdot = 0.
            results.capital_cost_model = Dummy()
            results.capital_cost_model.LCOE_brownfield = Dummy()
            results.capital_cost_model.LCOE_brownfield.LCOE = np.nan
        lcoe = results.capital_cost_model.LCOE_brownfield.LCOE * 1e6
        print(results.optMdot)
        print(lcoe)
        if np.isnan(lcoe):
            lcoe = 0.
        else:
            lcoe = lcoe[0]
        print(error_str)
        file.write(','.join([str(i) for i in [depth, permeability, results.optMdot, lcoe, """%s\n"""%error_str]]))
        # 1/0
file.close()
