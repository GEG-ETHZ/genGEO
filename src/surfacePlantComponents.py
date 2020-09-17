
from utils.fluidStates import FluidState

class DownHolePumpOutput(object):
    """DownHolePumpOutput."""
    pass

class DownHolePump(object):
    """DownHolePump."""

    def __init__(self, well, pump_depth, max_pump_dP):
        self.well = well
        self.pump_depth = pump_depth
        self.max_pump_dP = max_pump_dP

    def solve(self, initial_state, m_dot, time_years, P_inj_surface):

        # Now pumping and second well
        dP_pump = 0
        dP_surface = -100
        d_dP_pump = 1e5 # 1e5 is one bar increments to increase pressure by pumping

        while dP_surface <= -100:
            try:
                if dP_pump > self.max_pump_dP:
                    dP_pump = self.max_pump_dP

                P_prod_pump_out = initial_state.P_Pa + dP_pump
                T_prod_pump_out = initial_state.T_C
                temp_at_pump_depth = self.well.T_e_initial + self.well.dT_dz * self.pump_depth

                state = FluidState.getStateFromPT(P_prod_pump_out, T_prod_pump_out, self.well.fluid)
                state = self.well.solve(state, m_dot, time_years)

                dP_surface = (state.P_Pa - P_inj_surface)
                if dP_pump == self.max_pump_dP:
                    break

                if dP_surface < 0:
                    d_dP_pump = -1 * dP_surface
                    dP_pump = dP_pump + d_dP_pump
            except:
                # print(val)
                # 1/0
                raise ValueError('asd')
               #  # Only catch problems of flashing fluid
               # if (strcmp(ME.identifier,'SemiAnalytic:BelowSaturationPressure'))
               #      dP_pump = dP_pump + d_dP_pump
               # else:
               #      rethrow(ME)

        return state

    def gatherOutput(self):
        output = DownHolePumpOutput()
        return output
