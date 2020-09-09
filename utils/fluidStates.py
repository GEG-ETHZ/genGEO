from CoolProp.CoolProp import PropsSI
from utils.globalConstants import globalConstants

class FluidState(object):
    """FluidState pulls fluid states from coolprop and provides methods for
        - single properties
        - fluid states with all properties.
    Units are:
        - P [Pa]
        - T [C]
        - h [J/kg]
        - cp [J/K]
        - rho [kg/m^3]
        - S [J/K]
        """
    @staticmethod
    def getHFromPT(P_Pa, T_C, fluid):
        return PropsSI('HMASS', 'P', P_Pa, 'T', T_C + globalConstants.kelvin2celsius, fluid)

    @staticmethod
    def getRhoFromPT(P_Pa, T_C, fluid):
        return PropsSI('DMASS', 'P', P_Pa, 'T', T_C + globalConstants.kelvin2celsius, fluid)

    @staticmethod
    def getCpFromPT(P_Pa, T_C, fluid):
        return PropsSI('CPMASS', 'P', P_Pa, 'T', T_C + globalConstants.kelvin2celsius, fluid)

    @staticmethod
    def getSFromPT(P_Pa, T_C, fluid):
        return PropsSI('SMASS', 'P', P_Pa, 'T', T_C + globalConstants.kelvin2celsius, fluid)

    @staticmethod
    def getVFromPT(P_Pa, T_C, fluid):
        return PropsSI('V', 'P', P_Pa, 'T', T_C + globalConstants.kelvin2celsius, fluid)

    @staticmethod
    def getTFromPh(P_Pa, h_Jkg, fluid):
        return PropsSI('T', 'P', P_Pa, 'HMASS', h_Jkg, fluid) - globalConstants.kelvin2celsius

    @staticmethod
    def getRhoFromPh(P_Pa, h_Jkg, fluid):
        return PropsSI('DMASS', 'P', P_Pa, 'HMASS', h_Jkg, fluid)

    @staticmethod
    def getCpFromPh(P_Pa, h_Jkg, fluid):
        return PropsSI('CPMASS', 'P', P_Pa, 'HMASS', h_Jkg, fluid)

    @staticmethod
    def getSFromPh(P_Pa, h_Jkg, fluid):
        return PropsSI('SMASS', 'P', P_Pa, 'HMASS', h_Jkg, fluid)

    @staticmethod
    def getVFromPh(P_Pa, h_Jkg, fluid):
        return PropsSI('V', 'P', P_Pa, 'HMASS', h_Jkg, fluid)

    @staticmethod
    def getPFromTQ(T_C, Q, fluid):
        return PropsSI('P', 'T', T_C + globalConstants.kelvin2celsius, 'Q', Q, fluid)

    @staticmethod
    def getStateFromPT(P_Pa, T_C, fluid):
        st = FluidState()
        st.P_Pa = P_Pa
        st.T_C = T_C
        st.fluid = fluid
        st.h_Jkg = st.getHFromPT(P_Pa, T_C, fluid)
        st.rho_kgm3 = st.getRhoFromPT(P_Pa, T_C, fluid)
        st.cp_JK = st.getCpFromPT(P_Pa, T_C, fluid)
        st.v_Pas = st.getVFromPT(P_Pa, T_C, fluid)
        st.S_JK = st.getSFromPT(P_Pa, T_C, fluid)
        return st

    @staticmethod
    def getStateFromPh(P_Pa, h_Jkg, fluid):
        st = FluidState()
        st.P_Pa = P_Pa
        st.h_Jkg = h_Jkg
        st.fluid = fluid
        st.T_C = st.getTFromPh(P_Pa, h_Jkg, fluid)
        st.rho_kgm3 = st.getRhoFromPh(P_Pa, h_Jkg, fluid)
        st.cp_JK = st.getCpFromPh(P_Pa, h_Jkg, fluid)
        st.v_Pas = st.getVFromPh(P_Pa, h_Jkg, fluid)
        st.S_JK = st.getSFromPh(P_Pa, h_Jkg, fluid)
        return st


if __name__ == '__main__':
    import numpy as np
    P = np.array([1.e6, 1.e7])
    T = np.array([50., 60.])
    state = FluidState.getStateFromPT(P, T, 'water')
    print(state.h_Jkg)
    print(state.rho_kgm3)
    print(state.cp_JK)
    print(state.S_JK)
    state = FluidState.getStateFromPh(P, state.h_Jkg, 'water')
    print(state.T_C)
    print(state.rho_kgm3)
    print(state.cp_JK)
    print(state.S_JK)
    print(FluidState.getPFromTQ(T, 0, 'water'))
