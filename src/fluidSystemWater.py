import numpy as np

from utils.fluidStateFromPT import FluidStateFromPT
from utils.solver import Solver

class FluidSystemWaterOutput(object):
    """FluidSystemWaterOutput."""
    pass

class FluidSystemWater(object):
    """FluidSystemWater provides methods to compute the water fluid cycle."""

    def __init__(self, params):

        self.params = params

        self.fluid = 'Water'
        self.injection_well = None
        self.reservoir = None
        self.production_well1 = None
        self.pump = None
        self.pp = None

    def solve(self, initial_state):
        P_system_min = FluidStateFromPT.getPFromTQ(self.params.T_reservoir(), 0, self.fluid) + 1e5
        initial_state.P_Pa_in = initial_state.P_Pa() + P_system_min
        injection_state = FluidStateFromPT(initial_state.P_Pa(), initial_state.T_C(), initial_state.fluid)
        # Find necessary injection pressure
        dP_downhole = np.nan
        dP_solver = Solver()
        dP_loops = 1
        stop =  False
        while np.isnan(dP_downhole) or abs(dP_downhole) > 10e3:
            injection_well_state    = self.injection_well.solve(injection_state)
            reservoir_state         = self.reservoir.solve(injection_well_state.state)

            # if already at P_system_min, stop looping
            if stop:
                break

            # find downhole pressure difference (negative means overpressure)
            dP_downhole = self.params.P_reservoir() - reservoir_state.state.P_Pa()
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

        if reservoir_state.state.P_Pa() >= self.params.P_reservoir_max():
            raise Exception('GenGeo::FluidSystemWater:ExceedsMaxReservoirPressure - '
                        'Exceeds Max Reservoir Pressure of %.3f MPa!'%(self.params.P_reservoir_max()/1e6))

        production_well1_state  = self.production_well1.solve(reservoir_state.state)
        production_well2_state  = self.pump.solve(production_well1_state.state, injection_state.P_Pa())
        pp_state                = self.pp.solve(production_well2_state.well.state)
        return pp_state.state

    def gatherOutput(self):
        output = FluidSystemWaterOutput()
        output.injection_well       = self.injection_well.gatherOutput()
        output.reservoir            = self.reservoir.gatherOutput()
        output.production_well1     = self.production_well1.gatherOutput()
        output.pump                 = self.pump.gatherOutput()
        output.production_well2     = self.pump.well.gatherOutput()
        output.pp                   = self.pp.gatherOutput()
        return output
