
from utils.fluidStateProperties import FluidStateProperties

class FluidStateFromPT(FluidStateProperties):
    """FluidStateFromPT."""

    def __init__(self, P_Pa, T_C, fluid):
        super(FluidStateProperties, self)
        self.P_Pa_in = P_Pa
        self.T_C_in = T_C
        self.fluid = fluid

    def P_Pa(self):
        return self.P_Pa_in

    def T_C(self):
        return self.T_C_in


if __name__ == '__main__':
    import numpy as np
    P = np.array([1.e6, 1.e7])
    T = np.array([50., 60.])
    state = FluidStateFromPT(P, T, 'water')
    print(state.P_Pa())
    print(state.T_C())
    print(state.h_Jkg())
    print(state.rho_kgm3())
    print(state.cp_JK())
    print(state.S_JK())
    print(state.mu_Pas())
