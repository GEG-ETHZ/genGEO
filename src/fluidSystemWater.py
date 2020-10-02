
# from src.porousReservoir import PorousReservoir
# from tests.testAssertion import testAssert

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
        self.pp = None

    def solve(self, initial_state, m_dot, time_years):
        injection_state = initial_state
        # Find necessary injection pressure
        # Reset to 1 to make sure it is re-run at least once every time
        dP_downhole = 1
        while dP_downhole > 0:
            injection_well_state    = self.injection_well.solve(injection_state, m_dot, time_years)
            reservoir_state         = self.reservoir.solve(injection_well_state, m_dot, time_years)
            # find downhole pressure difference (negative means overpressure)
            dP_downhole = self.reservoir.P_reservoir - reservoir_state.P_Pa()
            if dP_downhole > 0.:
                injection_state.P_Pa_in = injection_state.P_Pa() + dP_downhole

        if reservoir_state.P_Pa() >= self.reservoir.P_reservoir_max:
            raise ValueError('FluidSystemWater:ExceedsMaxReservoirPressure - '
                        'Exceeds Max Reservoir Pressure of %.3f MPa!'%(self.reservoir.P_reservoir_max/1e6))

        production_well1_state  = self.production_well1.solve(reservoir_state, m_dot, time_years)
        production_well2_state  = self.pump.solve(production_well1_state, m_dot, time_years, injection_state.P_Pa())
        pp_state                = self.pp.solve(production_well2_state, m_dot, time_years)
        return pp_state

    def gatherOutput(self):
        output = FluidSystemWaterOutput()
        output.injection_well       = self.injection_well.gatherOutput()
        output.reservoir            = self.reservoir.gatherOutput()
        output.production_well1     = self.production_well1.gatherOutput()
        output.pump                 = self.pump.gatherOutput()
        output.production_well2     = self.pump.well.gatherOutput()
        output.pp                   = self.pp.gatherOutput()
        return output
