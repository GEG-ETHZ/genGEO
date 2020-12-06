
class SimulationParameters(object):
    """SimulationParameters provides physical properties of the system."""

    def __init__(self,
                working_fluid = None,
                orc_fluid = None,
                m_dot_IP = 1.,
                time_years = 1.,
                # subsurface model
                depth = 2500.,
                pump_depth = 500.,
                well_radius = 0.205,
                well_multiplier = 1.,
                well_spacing = 707.,
                well_length = None,
                monitoring_well_length = None,
                monitoring_well_radius = 0.108,
                dT_dz = 0.035,
                silica_precipitation = False,
                # T_e_initial = 15.,
                T_surface_rock = 15,
                T_ambient_C = 15.,
                reservoir_thickness = 100.,
                permeability = 1.0e-15 * 15000 / 100., # permeability = transmissivity / thickness
                reservoir_configuration = '5spot',
                has_surface_gathering_system = True,
                # power plant model
                max_pump_dP = 10.e6,
                eta_pump = 0.75,
                dT_approach = 7.,
                dT_pinch = 5.,
                eta_pump_orc = 0.9,
                eta_turbine_orc = 0.8,
                eta_pump_co2 = 0.9,
                eta_turbine_co2 = 0.78,
                cooling_mode = 'Wet',
                # cost model
                cost_year = 2019,
                success_rate = 0.95,
                well_cost_N = 1,
                pipe_cost_N = 0,
                F_OM = 0.045,
                discount_rate = 0.096,
                lifetime = 25,
                capacity_factor = 0.85,
                opt_mode = None,
                # physical properties
                g = 9.81,                       # m/s**2
                rho_rock = 2650.,               # kg/m**3
                c_rock = 1000.,                 # J/kg-K
                k_rock = 2.1,                   # W/m-K
                useWellboreHeatLoss = True,     # bool
                N_dx = 100,                     # nb well segments
                # Friction factor
                epsilon = 55 * 1e-6             # um
                    ):

        self.working_fluid = working_fluid
        self.orc_fluid = orc_fluid
        self.m_dot_IP = m_dot_IP
        self.time_years = time_years
        self.depth = depth
        self.pump_depth = pump_depth
        self.well_radius = well_radius
        self.well_multiplier = well_multiplier
        self.well_spacing = well_spacing
        self.monitoring_well_radius = monitoring_well_radius
        self.dT_dz = dT_dz
        self.silica_precipitation = silica_precipitation
        # self.T_e_initial = T_e_initial
        self.T_surface_rock = T_surface_rock
        self.T_ambient_C = T_ambient_C
        self.reservoir_thickness = reservoir_thickness
        self.permeability = permeability
        self.reservoir_configuration = reservoir_configuration
        self.has_surface_gathering_system = has_surface_gathering_system
        self.max_pump_dP = max_pump_dP
        self.eta_pump = eta_pump
        self.dT_approach = dT_approach
        self.dT_pinch = dT_pinch
        self.eta_pump_orc = eta_pump_orc
        self.eta_turbine_orc = eta_turbine_orc
        self.eta_pump_co2 = eta_pump_co2
        self.eta_turbine_co2 = eta_turbine_co2
        self.cooling_mode = cooling_mode
        if well_length:
            self.well_length = well_length
        else:
            self.well_length = self.depth
        if monitoring_well_length:
            self.monitoring_well_length = monitoring_well_length
        else:
            self.monitoring_well_length = self.depth
        self.cost_year = cost_year
        self.success_rate = success_rate
        self.well_cost_N = well_cost_N
        self.pipe_cost_N = pipe_cost_N
        self.F_OM = F_OM
        self.discount_rate = discount_rate
        self.lifetime = lifetime
        self.capacity_factor = capacity_factor
        self.opt_mode = opt_mode
        self.g = g
        self.rho_rock = rho_rock
        self.c_rock = c_rock
        self.k_rock = k_rock
        self.useWellboreHeatLoss = useWellboreHeatLoss
        self.N_dx = N_dx
        self.epsilon = epsilon
