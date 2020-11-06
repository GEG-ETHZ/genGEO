
from utils.fluidStates import FluidState

class FluidStateProperties(FluidState):
    """FluidStateProperties."""

    def __init__(self):
        super().__init__()

    def rho_kgm3(self):
        return self.getRhoFromPh(self.P_Pa(), self.h_Jkg(), self.fluid)

    def cp_JK(self):
        return self.getCpFromPh(self.P_Pa(), self.h_Jkg(), self.fluid)

    def mu_Pas(self):
        return self.getMuFromPh(self.P_Pa(), self.h_Jkg(), self.fluid)

    def S_JK(self):
        return self.getSFromPh(self.P_Pa(), self.h_Jkg(), self.fluid)

    def h_Jkg(self):
        return self.getHFromPT(self.P_Pa(), self.T_C(), self.fluid)

    def T_C(self):
        return self.getTFromPh(self.P_Pa(), self.h_Jkg(), self.fluid)
