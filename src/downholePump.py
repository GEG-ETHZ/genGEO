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
from models.simulationParameters import SimulationParameters
from utils.frictionFactor import frictionFactor

class DownHolePumpOutput(object):
    """DownHolePumpOutput."""
    w_pump = None

class DownHolePump(object):
    """DownHolePump."""

    def __init__(self, well, params = None, **kwargs):
        self.well = well
        self.params = params
        if self.params == None:
            self.params = SimulationParameters(**kwargs)
        self.computeSurfacePipeFrictionFactor()

    def computeSurfacePipeFrictionFactor(self):
        self.ff_m_dot = self.params.m_dot_IP
        initial_P = 1e6 +  self.params.P_system_min()
        initial_T = 60.
        initial_h = FluidState.getStateFromPT(initial_P, initial_T, self.params.working_fluid).h_Jkg
        self.friction_factor = frictionFactor(self.params.well_radius, initial_P, initial_h, \
                            self.params.m_dot_IP, self.params.working_fluid, self.params.epsilon)

    def solve(self, initial_state, P_inj_surface):

        if self.ff_m_dot != self.params.m_dot_IP:
            self.computeSurfacePipeFrictionFactor()

        # initilize object with well output
        results = DownHolePumpOutput()

        # Pumping and second well
        d_dP_pump = 1e5
        dP_pump = 0
        dP_surface = np.nan
        dP_loops = 1
        dP_Solver = Solver()

        while np.isnan(dP_surface) or abs(dP_surface) > 100:
            try:
                if dP_pump > self.params.max_pump_dP:
                    dP_pump = self.params.max_pump_dP

                P_prod_pump_out = initial_state.P_Pa + dP_pump
                T_prod_pump_out = initial_state.T_C
                temp_at_pump_depth = self.well.T_e_initial + self.params.dT_dz * self.params.pump_depth

                state_in = FluidState.getStateFromPT(P_prod_pump_out, T_prod_pump_out, self.params.working_fluid)
                results.well = self.well.solve(state_in)

                # calculate surface pipe friction loss
                if self.params.has_surface_gathering_system:
                    dP_surface_pipes = self.friction_factor * self.params.well_spacing / \
                    (2 * self.params.well_radius)**5 * 8 * self.params.m_dot_IP**2 / results.well.state.rho_kgm3 / np.pi**2
                else:
                    dP_surface_pipes = 0.

                dP_surface = (results.well.state.P_Pa - P_inj_surface - dP_surface_pipes)
                if dP_pump == self.params.max_pump_dP:
                    break

                # No pumping is needed in this system
                if dP_pump == 0 and dP_surface >= 0:
                    break

                dP_pump = dP_Solver.addDataAndEstimate(dP_pump, dP_surface)
                if np.isnan(dP_pump):
                    dP_pump = -1 * dP_surface

                # Pump can't be less than zero
                if dP_pump < 0:
                    dP_pump = 0

                # Warn against excessive loops
                if dP_loops > 10:
                    print('Warning::DownHolePump:dP_loops is large: %s'%dP_loops)
                dP_loops += 1

            except ValueError as error:
                # Only catch problems of flashing fluid
                if str(error).find('GenGeo::SemiAnalyticalWell:BelowSaturationPressure') > -1:
                    dP_pump = dP_pump + d_dP_pump
                else:
                   raise error
        # if pump pressure greater than allowable, throw error
        if dP_pump >= self.params.max_pump_dP:
            raise Exception('GenGeo::DownHolePump:ExceedsMaxProductionPumpPressure - '
            'Exceeds Max Pump Pressure of %.3f MPa!' %(self.params.max_pump_dP/1e6))

        results.w_pump = 0
        results.dP_surface_pipes = dP_surface_pipes
        if dP_pump > 0:
            results.w_pump = (initial_state.h_Jkg - state_in.h_Jkg) / self.params.eta_pump

        return results
