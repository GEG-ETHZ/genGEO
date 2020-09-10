from utils.fluidStates import FluidState

class ORCCycleResults(object):
    """ORCCycleResults."""
    def __init__(self):
        self.q_preheater = None
        self.q_boiler = None
        self.w_turbine = None
        self.q_desuperheater = None
        self.q_condenser = None
        self.w_pump = None
        self.w_cooler = None
        self.w_condenser = None
        self.w_net = None
        self.dT_range_CT = None
        self.dT_LMTD_preheater = None
        self.dT_LMTD_boiler = None
        self.end_T_C = None
        self.fluid =  'Water'

    def finalState(self):
        return FluidState.getStateFromTQ(self.end_T_C, 0, self.fluid)
