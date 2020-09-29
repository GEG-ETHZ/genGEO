
from src.energyConversion import EnergyConversion

class FullSystemOutput(object):
    """FullSystemOutput."""
    pass


class FullSystem(object):
    """FullSystem."""

    def __init__(self, fluid_system_solver, capital_cost_model):
        self.fluid_system_solver = fluid_system_solver
        self.capital_cost_model = capital_cost_model

    def solve(self, m_dot, time_years):
        self.m_dot = m_dot
        self.fluid_system_solver.solve(m_dot, time_years)
        self.energy_results = EnergyConversion.gatherOutput(self.m_dot, self.fluid_system_solver)
        self.capital_cost_model.energy_results = self.energy_results
        self.capital_cost_model.fluid_system = self.fluid_system_solver.fluid_system
        self.capital_cost_model.solve()

    def gatherOutput(self):
        output = FullSystemOutput()
        output.fluid_system_solver = self.fluid_system_solver.gatherOutput()
        output.energy_results = EnergyConversion.gatherOutput(self.m_dot, self.fluid_system_solver)
        output.capital_cost_model = self.capital_cost_model.gatherOutput()
        return output
