
from utils.fluidStateProperties import FluidStateProperties

class FluidStateFromPh(FluidStateProperties):
    """FluidStateFromPh."""

    def __init__(self, P_Pa, h_Jkg, fluid):
        super(FluidStateProperties, self).__init__()
        self.P_Pa_init = P_Pa
        self.h_Jkg_init = h_Jkg
        self.fluid = fluid

    def P_Pa(self):
        return self.P_Pa_init

    def h_Jkg(self):
        return self.h_Jkg_init


if __name__ == '__main__':
    import numpy as np
    from fluidStates import FluidState
    P = np.array([1.e6, 1.e7])
    T = np.array([50., 60.])
    h_Jkg = FluidState.getHFromPT(P, T, 'water')
    state = FluidStateFromPh(P, h_Jkg, 'water')
    print(state.P_Pa())
    print(state.T_C())
    print(state.h_Jkg())
    print(state.rho_kgm3())
    print(state.cp_JK())
    print(state.S_JK())
    print(state.v_Pas())
