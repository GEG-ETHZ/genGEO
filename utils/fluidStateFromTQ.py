
from utils.fluidStateProperties import FluidStateProperties

class FluidStateFromTQ(FluidStateProperties):
    """FluidStateFromTQ."""

    def __init__(self, T_C, Q, fluid):
        super(FluidStateProperties, self).__init__()
        self.Q_init = Q
        self.T_C_init = T_C
        self.fluid = fluid

    def P_Pa(self):
        return self.getPFromTQ(self.T_C(), self.Q_init, self.fluid)

    def h_Jkg(self):
        return self.getHFromTQ(self.T_C(), self.Q_init, self.fluid)

    def T_C(self):
        return self.T_C_init

if __name__ == '__main__':
    import numpy as np
    T = np.array([50., 60.])
    state = FluidStateFromTQ(T, 0, 'water')
    print(state.P_Pa())
    print(state.T_C())
    print(state.h_Jkg())
    print(state.rho_kgm3())
    print(state.cp_JK())
    print(state.S_JK())
    print(state.v_Pas())
