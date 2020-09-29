

class EnergyConversionORC(object):
    """EnergyConversionORC."""

    @staticmethod
    def gatherOutput(m_dot, input):
        N_IP_multiplier = input.fluid_system.injection_well.wellMultiplier
        input = input.gatherOutput()
        ec = EnergyConversionORC()
        ec.Q_preheater = m_dot * input.pp.q_preheater
        ec.Q_recuperator = m_dot * input.pp.q_recuperator
        ec.Q_boiler = m_dot * input.pp.q_boiler
        ec.Q_desuperheater = m_dot * input.pp.q_desuperheater
        ec.Q_condenser = m_dot * input.pp.q_condenser
        ec.W_turbine = m_dot * input.pp.w_turbine
        ec.W_pump = m_dot * input.pp.w_pump
        ec.W_cooler = m_dot * input.pp.w_cooler
        ec.W_condenser = m_dot * input.pp.w_condenser

        ec.Q_fluid = ec.Q_preheater + ec.Q_boiler

        ec.W_downhole_pump = m_dot * input.pump.w_pump
        ec.W_net = m_dot * input.pp.w_net + ec.W_downhole_pump

        ec.Q_preheater_total = N_IP_multiplier * ec.Q_preheater
        ec.Q_boiler_total = N_IP_multiplier * ec.Q_boiler
        ec.Q_recuperator_total = N_IP_multiplier * ec.Q_recuperator
        ec.Q_desuperheater_total = N_IP_multiplier * ec.Q_desuperheater
        ec.Q_condenser_total = N_IP_multiplier * ec.Q_condenser
        ec.W_turbine_total = N_IP_multiplier * ec.W_turbine
        ec.W_pump_orc_total = N_IP_multiplier * ec.W_pump
        ec.W_pump_prod_total = N_IP_multiplier * ec.W_downhole_pump
        ec.W_net_total = N_IP_multiplier * ec.W_net
        return ec

class EnergyConversionCPG(object):
    """EnergyConversionCPG."""

    @staticmethod
    def gatherOutput(m_dot, input):
        N_IP_multiplier = input.injection_well.wellMultiplier
        input = input.gatherOutput()
        ec = EnergyConversionCPG()
        ec.Q_preheater = m_dot * input.pp.q_preheater
        ec.Q_recuperator = m_dot * input.pp.q_recuperator
        ec.Q_boiler = m_dot * input.pp.q_boiler
        ec.Q_desuperheater = 0. # m_dot * input.pp.q_desuperheater
        ec.Q_condenser = m_dot * input.pp.q_condenser
        ec.W_turbine = m_dot * input.pp.w_turbine
        ec.W_pump = m_dot * input.pp.w_pump
        ec.W_cooler = m_dot * input.pp.w_cooler
        ec.W_condenser = m_dot * input.pp.w_condenser

        ec.W_net = m_dot * input.pp.w_net

        ec.Q_preheater_total = N_IP_multiplier * ec.Q_preheater
        ec.Q_boiler_total = N_IP_multiplier * ec.Q_boiler
        ec.Q_recuperator_total = N_IP_multiplier * ec.Q_recuperator
        ec.Q_desuperheater_total = N_IP_multiplier * ec.Q_desuperheater
        ec.Q_condenser_total = N_IP_multiplier * ec.Q_condenser
        ec.W_turbine_total = N_IP_multiplier * ec.W_turbine
        ec.W_pump_total = N_IP_multiplier * ec.W_pump
        ec.W_net_total = N_IP_multiplier * ec.W_net
        return ec
