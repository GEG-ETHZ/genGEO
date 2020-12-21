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
import copy
import numpy as np

class CapitalCostSystemOutput(object):
    """CapitalCostSystemOutput."""

class CapitalCostSystem(object):
    """CapitalCostSystem."""

    def __init__(self):
        self.energy_results = None
        self.fluid_system = None
        self.CapitalCost_SurfacePlant = None
        self.CapitalCost_SurfacePipe = None
        self.CapitalCost_Production_Well = None
        self.CapitalCost_Injection_Well = None
        self.CapitalCost_Wellfield = None
        self.CapitalCost_Exploration = None
        self.CapitalCost_Stimulation = None
        self.lcoe_model = None

    def solve(self):

        results = CapitalCostSystemOutput()

        lcoe_model_greenfield = copy.deepcopy(self.lcoe_model)
        lcoe_model_brownfield = copy.deepcopy(self.lcoe_model)

        results.C_surface_plant = self.CapitalCost_SurfacePlant.solve(self.energy_results, self.fluid_system)
        results.C_gathering_system = self.CapitalCost_SurfacePipe
        results.C_wells_production = self.CapitalCost_Production_Well
        results.C_wells_injection = self.CapitalCost_Injection_Well
        results.C_wellfield = self.CapitalCost_Wellfield
        results.C_exploration = self.CapitalCost_Exploration
        results.C_stimulation = self.CapitalCost_Stimulation

        C_surface_plant_plant = results.C_surface_plant.C_plant

        results.C_greenfield = np.sum([C_surface_plant_plant, results.C_gathering_system, results.C_wells_production, results.C_wells_injection, results.C_wellfield, results.C_exploration, results.C_stimulation])
        results.C_brownfield = np.sum([C_surface_plant_plant, results.C_gathering_system, results.C_wells_production])

        results.LCOE_greenfield = lcoe_model_greenfield.solve(results.C_greenfield, self.energy_results)
        results.LCOE_brownfield = lcoe_model_brownfield.solve(results.C_brownfield, self.energy_results)

        return results
