
from src.energyConversion import EnergyConversionORC
from src.fullSystem import FullSystemOutput

class FullSystemORC(object):
    """FullSystemORC."""

    def __init__(self, params, fluid_system_solver, capital_cost_model):
        self.params = params
        self.fluid_system_solver = fluid_system_solver
        self.capital_cost_model = capital_cost_model

    def solve(self):

        results = FullSystemOutput()

        results.fluid_system_solver = self.fluid_system_solver.solve()
        results.energy_results = EnergyConversionORC.compute(self.params, results.fluid_system_solver)
        self.capital_cost_model.energy_results = results.energy_results
        self.capital_cost_model.fluid_system = results.fluid_system_solver
        results.capital_cost_model = self.capital_cost_model.solve()

        return results
