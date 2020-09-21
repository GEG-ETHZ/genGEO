
# from src.porousReservoir import PorousReservoir

class FluidSystemWaterOutput(object):
    """FluidSystemWaterOutput."""
    pass

class FluidSystemWater(object):
    """FluidSystemWater provides methods to compute the water fluid cycle."""

    def __init__(self):
        self.fluid = 'Water'
        self.injection_well = None
        self.reservoir = None
        self.production_well1 = None
        self.pump = None
        self.oRC = None

    def solve(self, initial_state, m_dot, time_years):
        injection_state = initial_state

        # Find necessary injection pressure
        # Reset to 1 to make sure it is re-run at least once every time
        P_inj_surface = 1e6
        dP_downhole = 1
        while dP_downhole > 0:
            injection_well_state    = self.injection_well.solve(injection_state, m_dot, time_years)
            reservoir_state         = self.reservoir.solve(injection_well_state, m_dot, time_years)
            # find downhole pressure difference (negative means overpressure)
            dP_downhole = self.reservoir.P_reservoir - reservoir_state.P_Pa()
            if dP_downhole > 0.:
                P_inj_surface = P_inj_surface + dP_downhole

        production_well1_state  = self.production_well1.solve(reservoir_state, m_dot, time_years)
        production_well2_state  = self.pump.solve(production_well1_state, m_dot, time_years, P_inj_surface)
        oRC_state               = self.oRC.solve(injection_well_state, m_dot, time_years)

        return oRC_state

    def gatherOutput(self):
        output = FluidSystemWaterOutput()
        output.injection_well       = self.injection_well.gatherOutput()
        output.reservoir            = self.reservoir.gatherOutput()
        output.production_well1     = self.production_well1.gatherOutput()
        output.pump                 = self.pump.gatherOutput()
        output.production_well2     = self.pump.well.gatherOutput()
        output.oRC                  = self.oRC.gatherOutput()
        return output
