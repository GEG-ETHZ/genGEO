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

from CoolProp.CoolProp import PropsSI

def coolProps(outPropName, prop1Name, prop1Val, prop2Name, prop2Val, fluid):
    return PropsSI(outPropName, prop1Name, prop1Val, prop2Name, prop2Val, fluid)

P_crit_co2 = coolProps('PCRIT', "", 0, "", 0, 'CO2')
dP_tolerance = 5e3

def coolProp(outPropName, prop1Name, prop1Val, prop2Name, prop2Val, fluid):
    if fluid.lower() == 'co2':
        if prop1Name.lower() == 'p':
            P_lower = P_crit_co2 - dP_tolerance
            P_upper = P_crit_co2 + dP_tolerance
            if P_lower < prop1Val < P_upper:
                lower = coolProps(outPropName, 'P', P_lower, prop2Name, prop2Val, fluid)
                upper = coolProps(outPropName, 'P', P_upper, prop2Name, prop2Val, fluid)
                prop1Vals = np.array([P_lower, P_upper])
                prop2Vals = np.array([lower, upper])
                return np.interp(prop1Val, prop1Vals, prop2Vals)

    return coolProps(outPropName, prop1Name, prop1Val, prop2Name, prop2Val, fluid)
