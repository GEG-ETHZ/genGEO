
from src.energyConversion import EnergyConversionORC
from src.fullSystem import FullSystemOutput

class FullSystemORC(object):
    """FullSystemORC."""

    def __init__(self, params, fluid_system_solver, capital_cost_model):
        self.params = params
        self.fluid_system_solver = fluid_system_solver
        self.capital_cost_model = capital_cost_model

    def solve(self):
        self.fluid_system_solver.solve()
        self.energy_results = EnergyConversionORC.gatherOutput(self.params, self.fluid_system_solver.fluid_system)
        self.capital_cost_model.energy_results = self.energy_results
        self.capital_cost_model.fluid_system = self.fluid_system_solver.fluid_system
        self.capital_cost_model.solve()

    def gatherOutput(self):
        output = FullSystemOutput()
        output.fluid_system_solver = self.fluid_system_solver.gatherOutput()
        output.energy_results = self.energy_results
        output.capital_cost_model = self.capital_cost_model.gatherOutput()
        return output
