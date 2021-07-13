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
import os
import numpy as np

from src.fullSystemCPG import FullSystemCPG
from src.fullSystemSolver import FullSystemSolver

from models.simulationParameters import SimulationParameters

logTrans = np.arange(2., 8., 1.)
permeabilities = 1e-15 * 10. ** logTrans
depths = np.arange(1000, 8000, 1000)

# create output folder
output_folder = 'results'
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

output_file = open(os.path.join(output_folder, 'data_CO2.csv'), 'w')

# initialize parameters
params = SimulationParameters(working_fluid = 'co2', capacity_factor = 0.9)

# generate the full system
full_system = FullSystemCPG.getDefaultCPGSystem(params)
full_system_solver = FullSystemSolver(full_system)

# iterate over all depths and permeabilities and solve the system
for depth in depths:
    for permeability in permeabilities:
        print('Depth: ', depth)
        print('Permeability: ', permeability)
        params.depth = depth
        params.permeability = permeability / 100.

        try:
            output = full_system_solver.solve()
            lcoe_b = output.capital_cost_model.LCOE_brownfield.LCOE * 1e6
            lcoe_g = output.capital_cost_model.LCOE_greenfield.LCOE * 1e6
            power = output.energy_results.W_net / 1e6
            optMdot = output.optMdot
            error_str = ''

        except Exception as error:
            error_str = str(error).replace("\n", "").replace(",", " - ")
            lcoe_b = 0.
            lcoe_g = 0.
            power = 0.
            optMdot = 0.

        output_file.write(','.join([str(i) for i in [depth, permeability, optMdot, lcoe_b, lcoe_g, power, """%s\n"""%error_str]]))

output_file.close()
