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
from utils.coolPropInterface import coolProp
from utils.constantsAndPaths import ConversionConstants

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
    def __init__(self, prop1Name = None, prop1Val = None, prop2Name = None, prop2Val = None, fluid = None):
        self.prop1Name = prop1Name
        self.prop1Val = prop1Val
        self.prop2Name = prop2Name
        self.prop2Val = prop2Val
        self.fluid = fluid

    @property
    def P_Pa(self):
        return self.getCoolProp('P')

    @P_Pa.setter
    def P_Pa(self, P_Pa):
        self.setProp('P', P_Pa)

    @property
    def h_Jkg(self):
        return self.getCoolProp('HMASS')

    @h_Jkg.setter
    def h_Jkg(self, h_Jkg):
        self.setProp('HMASS', h_Jkg)

    @property
    def T_C(self):
        return self.getCoolProp('T') - ConversionConstants.kelvin2celsius

    @T_C.setter
    def T_C(self, T_C):
        self.setProp('T', T_C + ConversionConstants.kelvin2celsius)

    @property
    def s_JK(self):
        return self.getCoolProp('SMASS')

    @s_JK.setter
    def s_JK(self, s_JK):
        self.setProp('SMASS', s_JK)

    @property
    def rho_kgm3(self):
        return self.getCoolProp('DMASS')

    @property
    def cp_JK(self):
        return self.getCoolProp('CPMASS')

    @property
    def mu_Pas(self):
        return self.getCoolProp('V')

    def getCoolProp(self, prop):
        if self.prop1Name == prop:
            return self.prop1Val
        if self.prop2Name == prop:
            return self.prop2Val
        return coolProp(prop, self.prop1Name, self.prop1Val, self.prop2Name, self.prop2Val, self.fluid)

    def setProp(self, prop, value):
        if self.prop1Name == prop:
            self.prop1Val = value
        elif self.prop2Name == prop:
            self.prop2Val = value
        elif self.prop1Name == None:
            self.prop1Name = prop
            self.prop1Val = value
        elif self.prop2Name == None:
            self.prop2Name = prop
            self.prop2Val = value
        else:
            raise Exception('Property is over-specified.')

    @staticmethod
    def getStateFromPh(P_Pa, h_Jkg, fluid):
        return FluidState(prop1Name = 'P', prop1Val = P_Pa, prop2Name = 'HMASS', prop2Val = h_Jkg, fluid = fluid)

    @staticmethod
    def getStateFromPT(P_Pa, T_C, fluid):
        return FluidState(prop1Name = 'P', prop1Val = P_Pa, prop2Name = 'T', prop2Val = T_C + ConversionConstants.kelvin2celsius, fluid = fluid)

    @staticmethod
    def getStateFromPQ(P_Pa, Q, fluid):
        return FluidState(prop1Name = 'P', prop1Val = P_Pa, prop2Name = 'Q', prop2Val = Q, fluid = fluid)

    @staticmethod
    def getStateFromPS(P_Pa, s_JK, fluid):
        return FluidState(prop1Name = 'P', prop1Val = P_Pa, prop2Name = 'SMASS', prop2Val = s_JK, fluid = fluid)

    @staticmethod
    def getStateFromTQ(T_C, Q, fluid):
        return FluidState(prop1Name = 'T', prop1Val = T_C + ConversionConstants.kelvin2celsius, prop2Name = 'Q', prop2Val = Q, fluid = fluid)




    @staticmethod
    def getTcrit(fluid):
        return coolProp('TCRIT', "", 0, "", 0, fluid) - ConversionConstants.kelvin2celsius

    @staticmethod
    def getPcrit(fluid):
        return coolProp('PCRIT', "", 0, "", 0, fluid)
