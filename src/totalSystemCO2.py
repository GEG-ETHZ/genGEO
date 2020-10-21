
from src.semiAnalyticalWell import SemiAnalyticalWell
from src.porousReservoir import PorousReservoir
from src.subsurfaceComponents import DownHolePump
from src.oRCCycleTboil import ORCCycleTboil

from src.fluidSystemCO2 import FluidSystemCO2
from src.capitalCostSystem import CapitalCostSystem

from src.fullSystem import FullSystemCPG
from src.lCOESimple import LCOESimple
from src.capitalCostSurfacePipes import CapitalCostSurfacePipes
from src.capitalCostWell import CapitalCostWell
from src.capitalCostWellField import CapitalCostWellField
from src.capitalCostExploration import CapitalCostExploration
from src.capitalCostSurfacePlantCPG import CapitalCostSurfacePlantCPG
from src.fullSystemSolver import FullSystemSolverMinLCOEBrownfield

from utils.globalProperties import GlobalSimulationProperties
# from utils.globalConstants import globalConstants
# from utils.fluidStateFromPT import FluidStateFromPT

class TotalSystemCO2(object):
    """TotalSystemCO2."""

    def __init__(self,
                    working_fluid = 'CO2',
                    depth = 2500.,
                    well_radius = 0.205,
                    well_multiplier = 4.,
                    well_spacing = 707.,
                    monitoring_well_diameter = 0.216,
                    dT_dz = 0.035,
                    T_e_initial = 15.,
                    T_surface_rock = 15,
                    T_ambient_C = 15.,
                    reservoir_thickness = 100.,
                    permeability = 1.0e-15 * 15000 / 100., # permeability = transmissivity / thickness
                    reservoir_configuration = '5spot',
                    dT_approach = 7.,
                    eta_pump = 0.9,
                    eta_turbine = 0.78,
                    cooling_mode = 'Wet',
                    cost_year = 2019,
                    well_success_rate = 0.95,
                    F_OM = 0.045,
                    discount_rate = 0.096,
                    lifetime = 25,
                    capacity_factor = 0.85
                    ):


        # define global methods to be used in this tests
        gsp = GlobalSimulationProperties()

        fluid_system = FluidSystemCO2(T_ambient_C = T_ambient_C,
                                dT_approach = dT_approach,
                                eta_pump = eta_pump,
                                eta_turbine = eta_turbine,
                                coolingMode = cooling_mode)

        fluid_system.injection_well = SemiAnalyticalWell(params = gsp,
                                            dz_total = -depth,
                                            wellRadius = well_radius,
                                            wellMultiplier = well_multiplier,
                                            fluid = working_fluid,
                                            epsilon = 55 * 1e-6,
                                            dT_dz = dT_dz,
                                            T_e_initial = T_e_initial)

        fluid_system.reservoir = PorousReservoir(params = gsp,
                                        well_spacing = well_spacing,
                                        thickness = reservoir_thickness,
                                        permeability = permeability,
                                        T_surface_rock = T_surface_rock,
                                        depth = depth,
                                        dT_dz = dT_dz,
                                        wellRadius = well_radius,
                                        reservoirConfiguration = reservoir_configuration,
                                        fluid = working_fluid,
                                        modelPressureTransient = False,
                                        modelTemperatureDepletion = True)

        fluid_system.production_well = SemiAnalyticalWell(params = gsp,
                                            dz_total = depth,
                                            wellRadius = well_radius,
                                            wellMultiplier = well_multiplier,
                                            fluid = working_fluid,
                                            epsilon = 55 * 1e-6,
                                            dT_dz = dT_dz,
                                            T_e_initial = T_e_initial + dT_dz * depth)


        capital_cost_system = CapitalCostSystem(cost_year)
        capital_cost_system.CapitalCost_SurfacePlant = CapitalCostSurfacePlantCPG(cost_year)
        capital_cost_system.CapitalCost_SurfacePipe = CapitalCostSurfacePipes(N = 0)
        capital_cost_system.CapitalCost_Production_Well = CapitalCostWell.cO2Baseline(well_length = depth,
                                                                            well_diameter = well_radius * 2,
                                                                            success_rate = well_success_rate,
                                                                            cost_year = cost_year)
        capital_cost_system.CapitalCost_Injection_Well = CapitalCostWell.cO2Baseline(well_length = depth,
                                                                            well_diameter = well_radius * 2,
                                                                            success_rate = well_success_rate,
                                                                            cost_year = cost_year)
        capital_cost_system.CapitalCost_Wellfield = CapitalCostWellField.cO2MonitoringBaseline(N=1,
                                                                            monitoring_well_length = depth,
                                                                            monitoring_well_diameter = monitoring_well_diameter,
                                                                            cost_year=cost_year)
        capital_cost_system.CapitalCost_Exploration = CapitalCostExploration.cO2Baseline(N=1,
                                                                            well_length = depth,
                                                                            well_diameter = well_radius * 2,
                                                                            success_rate = well_success_rate,
                                                                            cost_year = cost_year)
        capital_cost_system.lcoe_model = LCOESimple(F_OM = F_OM,
                                                    discountRate = discount_rate,
                                                    Lifetime = lifetime,
                                                    CapacityFactor = capacity_factor)

        self.full_system = FullSystemCPG(fluid_system, capital_cost_system)
        self.full_system_solver = FullSystemSolver(self.full_system)

    def minimizeLCOEBrownfield(self, time_years = 1.):
        optMdot = self.full_system_solver.minimizeLCOEBrownfield(time_years = time_years)
        if optMdot:
            output = self.full_system.gatherOutput()
            output.optMdot = optMdot
            return output
        else:
            raise ValueError('No minimum found!')


    def minimizeLCOEGreenfield(self, time_years = 1.):
        optMdot = self.full_system_solver.minimizeLCOEGreenfield(time_years = time_years)
        if optMdot:
            output = self.full_system.gatherOutput()
            output.optMdot = optMdot
            return output
        else:
            raise ValueError('No minimum found!')


    def maximizePower(self, time_years = 1.):
        optMdot = self.full_system_solver.maximizePower(time_years = time_years)
        output = self.full_system.gatherOutput()
        output.optMdot = optMdot
        return output

    def solve(self, m_dot, time_years):
        self.full_system.solve(m_dot = m_dot, time_years = time_years)
        output = self.full_system.gatherOutput()
        return output
