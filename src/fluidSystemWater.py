import numpy as np

from utils.fluidStateFromPT import FluidStateFromPT
from utils.solver import Solver

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
        P_system_min = FluidStateFromPT.getPFromTQ(self.reservoir.reservoirT, 0, self.fluid) + 1e5
        initial_state.P_Pa_in = initial_state.P_Pa() + P_system_min
        injection_state = FluidStateFromPT(initial_state.P_Pa(), initial_state.T_C(), initial_state.fluid)
        # Find necessary injection pressure
        dP_downhole = np.nan
        dP_solver = Solver()
        dP_loops = 1
        stop =  False
        while np.isnan(dP_downhole) or abs(dP_downhole) > 10e3:
            injection_well_state    = self.injection_well.solve(injection_state, m_dot, time_years)
            reservoir_state         = self.reservoir.solve(injection_well_state, m_dot)

            # if already at P_system_min, stop looping
            if stop:
                break

            # find downhole pressure difference (negative means overpressure)
            dP_downhole = self.reservoir.P_reservoir - reservoir_state.P_Pa()
            injection_state.P_Pa_in = dP_solver.addDataAndEstimate(injection_state.P_Pa(), dP_downhole)

            if np.isnan(injection_state.P_Pa()):
                injection_state.P_Pa_in = initial_state.P_Pa() + dP_downhole

            if dP_loops > 10:
                print('Warning::FluidSystemWater:dP_loops is large: %s'%dP_loops)
            dP_loops += 1

            # Set Limits
            if injection_state.P_Pa() < P_system_min:
                # can't be below this temp or fluid will flash
                injection_state.P_Pa_in = P_system_min
                # switch stop to run injection well and reservoir once more
                stop = True

        if reservoir_state.P_Pa() >= self.reservoir.P_reservoir_max:
            raise Exception('FluidSystemWater:ExceedsMaxReservoirPressure - '
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
