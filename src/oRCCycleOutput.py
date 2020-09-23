
import numpy as np

class ORCCycleOutput(object):
    """ORCCycleOutput."""
    q_preheater = np.nan
    q_recuperator = np.nan
    q_boiler = np.nan
    q_desuperheater = np.nan
    q_condenser = np.nan
    w_turbine = np.nan
    w_pump = np.nan
    w_cooler = np.nan
    w_condenser = np.nan
    w_net = np.nan
    dT_range_CT = np.nan
    dT_LMTD_preheater = np.nan
    dT_LMTD_recuperator = np.nan
    dT_LMTD_boiler = np.nan
    state_out = np.nan
