
# from src.porousReservoir import PorousReservoir
from tests.testAssertion import testAssert

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
        self.orc = None

    def solve(self, initial_state, m_dot, time_years):
        injection_state = initial_state
        # Find necessary injection pressure
        # Reset to 1 to make sure it is re-run at least once every time
        dP_downhole = 1
        # print('#############   new  ############')
        while dP_downhole > 0:
            injection_well_state    = self.injection_well.solve(injection_state, m_dot, time_years)
            # print(*testAssert(injection_well_state.P_Pa(), 2.528328896329052e+07, 'test1_pressure'))
            # print(*testAssert(injection_well_state.T_C(), 62.110938581107170, 'test1_temp'))
            # breakpoint()
            reservoir_state         = self.reservoir.solve(injection_well_state, m_dot, time_years)
            # print(*testAssert(reservoir_state.P_Pa(), 2.509002066687885e+07, 'test2_pressure'))
            # print(*testAssert(reservoir_state.T_C(), 1.025000000000000e+02, 'test2_temp'))
            # breakpoint()
            # find downhole pressure difference (negative means overpressure)
            dP_downhole = self.reservoir.P_reservoir - reservoir_state.P_Pa()
            if dP_downhole > 0.:
                injection_state.P_Pa_in = injection_state.P_Pa() + dP_downhole
            # print(*testAssert(reservoir_state.P_Pa(), 2.509002066687885e+07, 'test3_pressure'))
            # print(*testAssert(reservoir_state.T_C(), 1.025000000000000e+02, 'test3_temp'))
            # breakpoint()

        if reservoir_state.P_Pa() >= self.reservoir.P_reservoir_max:
            raise ValueError('TotalAnalyticSystemWater:ExceedsMaxReservoirPressure - '
                                'Exceeds Max Reservoir Pressure of %.3f MPa!'%(self.reservoir.P_reservoir_max/1e6))

        production_well1_state  = self.production_well1.solve(reservoir_state, m_dot, time_years)
        # print(*testAssert(production_well1_state.P_Pa(), 6.111714723757768e+06, 'test4_pressure'))
        # print(*testAssert(production_well1_state.T_C(), 88.164719464305890, 'test4_temp'))
        # breakpoint()
        production_well2_state  = self.pump.solve(production_well1_state, m_dot, time_years, injection_state.P_Pa())
        # print(*testAssert(production_well2_state.P_Pa(), 1.352117133629285e+06, 'test5_pressure'))
        # print(*testAssert(production_well2_state.T_C(), 81.255757432571240, 'test5_temp'))
        # breakpoint()
        orc_state               = self.orc.solve(production_well2_state, m_dot, time_years)
        # print(*testAssert(orc_state.P_Pa(), 1.352117133629285e+06, 'test6_pressure'))
        # print(*testAssert(orc_state.T_C(), 54.038702389379260, 'test6_temp'))
        # breakpoint()

        return orc_state

    def gatherOutput(self):
        output = FluidSystemWaterOutput()
        output.injection_well       = self.injection_well.gatherOutput()
        output.reservoir            = self.reservoir.gatherOutput()
        output.production_well1     = self.production_well1.gatherOutput()
        output.pump                 = self.pump.gatherOutput()
        output.production_well2     = self.pump.well.gatherOutput()
        output.orc                  = self.oRC.gatherOutput()
        return output
