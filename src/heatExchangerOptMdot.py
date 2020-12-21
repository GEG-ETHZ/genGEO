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

from src.heatExchanger import heatExchanger

def heatExchangerOptMdot(T_1_in, P_1, fluid_1, T_2_in, P_2, fluid_2, dT_pinch, T_1_out_min, maximizeHeatFromStream = '2'):

    mdot_ratio = 1
    d_mdot_ratio = 0.2
    q_exchanged_old = np.nan
    peaks = 0

    while peaks < 4.:
        mdot_ratio = mdot_ratio + d_mdot_ratio
        # mdot_ratio = m_dot_2 / m_dot_1
        m_dot_1 = 1
        m_dot_2 = m_dot_1 * mdot_ratio

        # [Q_exchanged, dT_LMTD, T_1_out, T_2_out, T_1, T_2, Q] =

        heat_exchanger_results = heatExchanger(T_1_in, P_1, m_dot_1, fluid_1, T_2_in, P_2, m_dot_2, fluid_2, dT_pinch)
        heat_exchanger_results.q_exchanged_1 = heat_exchanger_results.Q_exchanged / m_dot_1
        heat_exchanger_results.q_exchanged_2 = heat_exchanger_results.Q_exchanged / m_dot_2

        if maximizeHeatFromStream == '1':
            q_exchanged = heat_exchanger_results.q_exchanged_1
        elif maximizeHeatFromStream == '2':
            q_exchanged = heat_exchanger_results.q_exchanged_2
        else:
            raise Exception('GenGeo::heatExchangerOptMdot:unknownOptimizeMethod - Unknown Optimize Method')

        # if temp too low, make sure mdot_ratio is increasing
        if not np.isnan(T_1_out_min) and heat_exchanger_results.T_1_out < T_1_out_min and d_mdot_ratio < 0:
            d_mdot_ratio = -0.21 * d_mdot_ratio
            peaks = peaks + 1
            q_exchanged = np.nan
        elif not np.isnan(q_exchanged) and not np.isnan(q_exchanged_old) and q_exchanged < q_exchanged_old:
            d_mdot_ratio = -0.21 * d_mdot_ratio
            peaks = peaks + 1

        q_exchanged_old = q_exchanged

    heat_exchanger_results.mdot_ratio = mdot_ratio
    heat_exchanger_results.m_dot_1 = m_dot_1
    heat_exchanger_results.m_dot_2 = m_dot_2


    return heat_exchanger_results
