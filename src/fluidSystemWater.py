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
import numpy as np

from utils.fluidState import FluidState
from utils.solver import Solver

from utils.frictionFactor import frictionFactor

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

        results = FluidSystemWaterOutput()

        injection_state = FluidState.getStateFromPT(initial_state.P_Pa, initial_state.T_C, initial_state.fluid)
        # Find necessary injection pressure
        dP_downhole = np.nan
        dP_solver = Solver()
        dP_loops = 1
        stop =  False

        while np.isnan(dP_downhole) or abs(dP_downhole) > 10e3:
            results.injection_well    = self.injection_well.solve(injection_state)
            results.reservoir         = self.reservoir.solve(results.injection_well.state)

            # if already at P_system_min, stop looping
            if stop:
                break

            # find downhole pressure difference (negative means overpressure)
            dP_downhole = self.params.P_reservoir() - results.reservoir.state.P_Pa
            injection_state.P_Pa = dP_solver.addDataAndEstimate(injection_state.P_Pa, dP_downhole)

            if np.isnan(injection_state.P_Pa):
                injection_state.P_Pa = initial_state.P_Pa + dP_downhole

            if dP_loops > 10:
                print('GenGeo::Warning:FluidSystemWater:dP_loops is large: %s'%dP_loops)
            dP_loops += 1

            # Set Limits
            if injection_state.P_Pa < self.params.P_system_min():
                # can't be below this temp or fluid will flash
                injection_state.P_Pa = self.params.P_system_min()
                # switch stop to run injection well and reservoir once more
                stop = True

        if results.reservoir.state.P_Pa >= self.params.P_reservoir_max():
            raise Exception('GenGeo::FluidSystemWater:ExceedsMaxReservoirPressure - '
                        'Exceeds Max Reservoir Pressure of %.3f MPa!'%(self.params.P_reservoir_max()/1e6))

        # Production Well (Lower, to production pump)
        results.production_well1  = self.production_well1.solve(results.reservoir.state)
        # Upper half of production well
        results.pump              = self.pump.solve(results.production_well1.state, injection_state.P_Pa)

        # Subtract surface frictional losses between production wellhead and surface plant
        ff = frictionFactor(self.params.well_radius, results.pump.well.state.P_Pa, results.pump.well.state.h_Jkg,
            self.params.m_dot_IP, self.params.working_fluid, self.params.epsilon)
        if self.params.has_surface_gathering_system == True:
            dP_surfacePipes = ff * self.params.well_spacing / (self.params.well_radius*2)**5 * 8 * self.params.m_dot_IP**2 / results.pump.well.state.rho_kgm3 / np.pi**2
        else:
            dP_surfacePipes = 0
        
        results.surface_plant_inlet = FluidState.getStateFromPh(
            results.pump.well.state.P_Pa - dP_surfacePipes,
            results.pump.well.state.h_Jkg,
            self.params.working_fluid)

        results.pp = self.pp.solve(results.surface_plant_inlet)
        return results
