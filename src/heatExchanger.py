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

from utils.fluidState import FluidState

class HeatExchangerResults(object):
    """HeatExchangerResults contains all results information of a heat exchanger."""

    def __init__(self):
        self.Q_exchanged =  0.
        self.dT_LMTD =  np.nan
        self.T_1_out =  0.
        self.T_2_out =  0.
        self.T_1 =  np.zeros(0)
        self.T_2 =  np.zeros(0)
        self.Q =  np.zeros(0)

def heatExchanger(T_1_in, P_1, m_dot_1, fluid_1, T_2_in, P_2, m_dot_2, fluid_2, dT_pinch):

    results = HeatExchangerResults()

    # Only tested where Temp1 < Temp2
    if T_1_in > T_2_in:
        # throw(MException('HeatExchanger:BadTemps','Temp 1 is greater than Temp 2!'));
        # No heat exchanged
        results.T_1_out = T_1_in
        results.T_2_out = T_2_in
        return results

    if m_dot_1 <= 0 or m_dot_2 <= 0:
        raise Exception('GenGeo::HeatExchanger:NegativeMassFlow - Negative Massflow in Heat Exchanger')

    increments = 20
    direction = np.sign(T_1_in - T_2_in)

    # Check the phase on incoming fluid
    P_crit_1 = FluidState.getPcrit(fluid_1)
    if P_1 < P_crit_1:
        T_sat_1 = FluidState.getStateFromPQ(P_1, 1, fluid_1).T_C
        if T_sat_1 == T_1_in or T_sat_1 == T_2_in:
            raise Exception('GenGeo::HeatExchanger:TwoPhaseFluid - Fluid 1 enters or leaves two-phase!')

    P_crit_2 = FluidState.getPcrit(fluid_2)
    if P_2 < P_crit_2:
        T_sat_2 = FluidState.getStateFromPQ(P_2, 1, fluid_2).T_C
        if T_sat_2 == T_1_in or T_sat_2 == T_2_in:
            raise Exception('GenGeo::HeatExchanger:TwoPhaseFluid - Fluid 2 enters or leaves two-phase!')

    h_1_in = FluidState.getStateFromPT(P_1, T_1_in, fluid_1).h_Jkg
    T_1_max = T_2_in
    h_1_max = FluidState.getStateFromPT(P_1, T_1_max, fluid_1).h_Jkg
    T_1_max_practical = T_2_in + direction * dT_pinch
    h_1_max_practical = FluidState.getStateFromPT(P_1, T_1_max_practical, fluid_1).h_Jkg
    h_2_in = FluidState.getStateFromPT(P_2, T_2_in, fluid_2).h_Jkg
    T_2_max = T_1_in;
    h_2_max = FluidState.getStateFromPT(P_2, T_2_max, fluid_2).h_Jkg
    T_2_max_practical = T_1_in - direction*dT_pinch;
    h_2_max_practical = FluidState.getStateFromPT(P_2, T_2_max_practical, fluid_2).h_Jkg

    Q_1_max = abs( m_dot_1 * (h_1_in - h_1_max) )
    Q_2_max = abs( m_dot_2 * (h_2_in - h_2_max) )
    Q_1_max_practical = abs( m_dot_1 * (h_1_in - h_1_max_practical) )
    Q_2_max_practical = abs( m_dot_2 * (h_2_in - h_2_max_practical) )

    if abs(Q_1_max) < abs(Q_2_max):
        # limitingFluid = 1;
        Q_max = Q_1_max
        Q_max_practical = Q_1_max_practical
    else:
        # limtingFluid = 2;
        Q_max = Q_2_max
        Q_max_practical = Q_2_max_practical

    results.Q_exchanged = Q_max_practical
    ddT_pinch = 1

    length      = increments + 1
    results.Q   = np.zeros(length)
    h_1         = np.zeros(length)
    results.T_1 = np.zeros(length)
    h_2         = np.zeros(length)
    results.T_2 = np.zeros(length)
    dT          = np.zeros(length)
    UA          = np.zeros(length)

    while ddT_pinch > 0.1:

        dQ = results.Q_exchanged / increments

        results.Q[0] = 0.
        h_1[0] = h_1_in
        results.T_1[0] = T_1_in
        h_2[0] = (Q_2_max - results.Q_exchanged) / m_dot_2 + h_2_max
        results.T_2[0] = FluidState.getStateFromPh(P_2, h_2[0], fluid_2).T_C
        dT[0] = direction * (results.T_1[0] - results.T_2[0])
        UA[0] = dQ / dT[0]

        for i in range (1, increments+1):
            results.Q[i] = results.Q[i-1] + dQ
            h_1[i] = h_1[i-1] + dQ / m_dot_1
            h_2[i] = h_2[i-1] + dQ / m_dot_2
            results.T_1[i] = FluidState.getStateFromPh(P_1, h_1[i], fluid_1).T_C
            results.T_2[i] = FluidState.getStateFromPh(P_2, h_2[i], fluid_2).T_C
            dT[i] = direction * (results.T_1[i] - results.T_2[i])
            UA[i] = dQ / dT[i]

        min_dT = min(dT)
        ddT_pinch = dT_pinch - min_dT

        # Adjust Q_exchanged
        # Use proportional error approach
        change = ddT_pinch / (T_2_in - T_1_in)
        results.Q_exchanged = (1-change) * results.Q_exchanged

    results.dT_LMTD = results.Q_exchanged / sum(UA)
    effectiveness = results.Q_exchanged / Q_max
    results.T_1_out = results.T_1[-1]
    results.T_2_out = results.T_2[0]

    # Check the phase on leaving fluid
    if P_1 < P_crit_1:
        T_sat_1 = FluidState.getStateFromPQ(P_1, 1, fluid_1).T_C
        if T_sat_1 > T_1_in and T_sat_1 < results.T_1_out:
            print('Caution: Fluid 1 is phase changing in heat exchanger')
        if T_sat_1 == results.T_1_out:
            raise Exception('GenGeo::HeatExchanger:TwoPhaseFluid - Fluid 1 leaves two-phase!')

    if P_2 < P_crit_2:
        T_sat_2 = FluidState.getStateFromPQ(P_2, 1, fluid_2).T_C
        if T_sat_2 > T_2_in and T_sat_2 < results.T_2_out:
            print('Caution: Fluid 2 is phase changing in heat exchanger')
        if T_sat_2 == results.T_2_out:
            raise Exception('GenGeo::HeatExchanger:TwoPhaseFluid - Fluid 2 leaves two-phase!')

    return results
