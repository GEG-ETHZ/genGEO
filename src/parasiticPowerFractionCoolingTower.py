from utils.constantsAndPaths import ConversionConstants


def computeProcess(a, b, c, d, T_ambient_C, dT_approach_CT, dT_range_CT):
    A = a * (1/dT_approach_CT)
    B = b * (T_ambient_C+ConversionConstants.kelvin2celsius)
    C = c * (T_ambient_C+ConversionConstants.kelvin2celsius)/dT_approach_CT
    D = d * (1/(dT_approach_CT + dT_range_CT))
    return A + B + C + D

def parasiticPowerFractionCoolingTower(T_ambient_C, dT_approach_CT, dT_range_CT, coolingMode):
    def processWet(process):
        if process == 'cooling':
            return computeProcess(1.2, 0., -3.79e-3, 1.95e-2, T_ambient_C, dT_approach_CT, dT_range_CT)
        elif process == 'condensing':
            return computeProcess(1.65, -6.24e-6, -5.03e-3, 0., T_ambient_C, dT_approach_CT, dT_range_CT)
        else:
            raise ValueError('parasiticPowerFraction:UnknownProcess - Unknown Process - use "cooling" or "condensing"')

    def processDry(process):
        if process == 'cooling':
            return computeProcess(1.2, 0., 0., 0., T_ambient_C, dT_approach_CT, dT_range_CT)
        elif process == 'condensing':
            return computeProcess(0.619, 0., 0., 0., T_ambient_C, dT_approach_CT, dT_range_CT)
        else:
            raise ValueError('parasiticPowerFraction:UnknownProcess - Unknown Process - use "cooling" or "condensing"')

    if coolingMode == 'Wet':
        return processWet
    elif coolingMode == 'Dry':
        return processDry
    else:
        raise ValueError('parasiticPowerFraction:UnknownCoolingMode - Unknown Cooling Mode')
