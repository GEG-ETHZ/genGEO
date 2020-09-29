import copy

class CapitalCostSystemWaterOutput(object):
    """CapitalCostSystemWaterOutput."""

class CapitalCostSystemWater(object):
    """CapitalCostSystemWater."""

    def __init__(self, cost_year):
        self.cost_year = cost_year
        self.energy_results = None
        self.fluid_system = None
        self.CapitalCost_SurfacePlant = None
        self.CapitalCost_SurfacePipe = None
        self.CapitalCost_Production_Well = None
        self.CapitalCost_Injection_Well = None
        self.CapitalCost_Wellfield = None
        self.CapitalCost_Exploration = None
        self.lcoe_model = None

    def solve(self):

        self.lcoe_model_greenfield = copy.deepcopy(self.lcoe_model)
        self.lcoe_model_brownfield = copy.deepcopy(self.lcoe_model)

        C_surfacePlant = self.CapitalCost_SurfacePlant.solve(self.energy_results, self.fluid_system).C_plant
        C_gatheringSystem = self.CapitalCost_SurfacePipe.solve(self.cost_year)
        C_wells_production = self.CapitalCost_Production_Well
        C_wells_injection = self.CapitalCost_Injection_Well
        C_wellfield = self.CapitalCost_Wellfield
        C_exploration = self.CapitalCost_Exploration
        C_stimulation = 0.

        self.capital_cost_greenfield = C_surfacePlant + C_gatheringSystem + C_wells_production + C_wells_injection + C_wellfield + C_exploration + C_stimulation
        self.capital_cost_brownfield = C_surfacePlant + C_gatheringSystem + C_wells_production

        self.lcoe_model_greenfield.solve(self.capital_cost_greenfield, self.energy_results)
        self.lcoe_model_brownfield.solve(self.capital_cost_brownfield, self.energy_results)

    def gatherOutput(self):
        output = CapitalCostSystemWaterOutput()
        output.C_greenfield = self.capital_cost_greenfield
        output.C_brownfield = self.capital_cost_brownfield
        output.LCOE_greenfield = self.lcoe_model_greenfield.gatherOutput()
        output.LCOE_brownfield = self.lcoe_model_brownfield.gatherOutput()
        return output
