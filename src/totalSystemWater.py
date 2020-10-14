
from src.semiAnalyticalWell import SemiAnalyticalWell
from src.porousReservoir import PorousReservoir
from src.subsurfaceComponents import DownHolePump
from src.oRCCycleTboil import ORCCycleTboil

from src.fluidSystemWater import FluidSystemWater
from src.fluidSystemWaterSolver import FluidSystemWaterSolver
from src.capitalCostSystem import CapitalCostSystem

from src.fullSystem import FullSystemORC
from src.lCOESimple import LCOESimple
from src.capitalCostSurfacePipes import CapitalCostSurfacePipes
from src.capitalCostWell import CapitalCostWell
from src.capitalCostWellField import CapitalCostWellField
from src.capitalCostExploration import CapitalCostExploration
from src.capitalCostSurfacePlantORC import CapitalCostSurfacePlantORC
from src.fullSystemSolver import FullSystemSolver

from utils.globalProperties import GlobalSimulationProperties

class TotalSystemWater(object):
    """TotalSystemWater."""

    def __init__(self,
                    working_fluid = 'water',
                    depth = 2500.,
                    pump_depth = 500.,
                    well_radius = 0.205,
                    well_multiplier = 4.,
                    well_spacing = 707.,
                    dT_dz = 0.035,
                    T_e_initial = 15.,
                    T_surface_rock = 15,
                    T_ambient_C = 15.,
                    reservoir_thickness = 100.,
                    permeability = 1.0e-15 * 15000 / 100., # permeability = transmissivity / thickness
                    reservoir_configuration = '5spot',
                    max_pump_dP = 10.e6,
                    eta_pump = 0.75,
                    dT_approach = 7.,
                    dT_pinch = 5.,
                    eta_pump_orc = 0.9,
                    eta_turbine_orc = 0.8,
                    cooling_mode = 'Wet',
                    orc_fluid = 'R600a',
                    cost_year = 2019,
                    well_success_rate = 0.95,
                    F_OM = 0.045,
                    discount_rate = 0.096,
                    lifetime = 25,
                    capacity_factor = 0.85
                    ):


        # define global methods to be used in this tests
        gsp = GlobalSimulationProperties()
        # alternative simulation properties for porous reservoir
        gsp2 = GlobalSimulationProperties()
        gsp2.k_rock = 2.1        #W/m/C
        gsp2.rho_rock = 2300     #kg/m^3
        gsp2.c_rock = 920.       #J/kg/C

        fluid_system = FluidSystemWater()
        fluid_system.injection_well = SemiAnalyticalWell(params = gsp,
                                            dz_total = -1 * depth,
                                            wellRadius = well_radius,
                                            wellMultiplier = well_multiplier,
                                            fluid = working_fluid,
                                            epsilon = 55 * 1e-6,
                                            dT_dz = dT_dz,
                                            T_e_initial = T_e_initial)

        fluid_system.reservoir = PorousReservoir(params = gsp2,
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

        fluid_system.production_well1 = SemiAnalyticalWell(params = gsp,
                                            dz_total = depth - pump_depth,
                                            wellRadius = well_radius,
                                            wellMultiplier = well_multiplier,
                                            fluid = working_fluid,
                                            epsilon = 55 * 1e-6,
                                            dT_dz = dT_dz,
                                            T_e_initial = T_e_initial + dT_dz * depth)

        prod_well2 = SemiAnalyticalWell(params = gsp,
                                            dz_total = pump_depth,
                                            wellRadius = well_radius,
                                            wellMultiplier = well_multiplier,
                                            fluid = working_fluid,
                                            epsilon = 55 * 1e-6,
                                            dT_dz = dT_dz,
                                            T_e_initial = T_e_initial + dT_dz * pump_depth)

        fluid_system.pump = DownHolePump(well = prod_well2,
                            pump_depth = pump_depth,
                            max_pump_dP = max_pump_dP,
                            eta_pump = eta_pump)

        fluid_system.pp = ORCCycleTboil(T_ambient_C = T_ambient_C,
                                dT_approach = dT_approach,
                                dT_pinch = dT_pinch,
                                eta_pump = eta_pump_orc,
                                eta_turbine = eta_turbine_orc,
                                coolingMode = cooling_mode,
                                orcFluid = orc_fluid)

        capital_cost_system = CapitalCostSystem(cost_year)
        capital_cost_system.CapitalCost_SurfacePlant = CapitalCostSurfacePlantORC(cost_year)
        capital_cost_system.CapitalCost_SurfacePipe = CapitalCostSurfacePipes(N = 0)
        capital_cost_system.CapitalCost_Production_Well = CapitalCostWell.waterBaseline(well_length = depth,
                                                                            well_diameter = well_radius * 2,
                                                                            success_rate = well_success_rate,
                                                                            cost_year = cost_year)
        capital_cost_system.CapitalCost_Injection_Well = CapitalCostWell.waterBaseline(well_length = depth,
                                                                            well_diameter = well_radius * 2,
                                                                            success_rate = well_success_rate,
                                                                            cost_year = cost_year)
        capital_cost_system.CapitalCost_Wellfield = CapitalCostWellField.water(cost_year)
        capital_cost_system.CapitalCost_Exploration = CapitalCostExploration.waterBaseline(well_length = depth,
                                                                            well_diameter = well_radius * 2,
                                                                            success_rate = well_success_rate,
                                                                            cost_year = cost_year)
        capital_cost_system.lcoe_model = LCOESimple(F_OM = F_OM,
                                                    discountRate = discount_rate,
                                                    Lifetime = lifetime,
                                                    CapacityFactor = capacity_factor)

        solver = FluidSystemWaterSolver(fluid_system)

        self.full_system = FullSystemORC(solver, capital_cost_system)

        self.full_system_solver = FullSystemSolver(self.full_system)

    def minimizeLCOEBrownfield(self, time_years = 1.):
        optMdot = self.full_system_solver.minimizeLCOEBrownfield(time_years = time_years)
        if optMdot:
            output = self.full_system.gatherOutput()
            output.optMdot = optMdot
            return output
        else:
            raise Exception('No minimum found!')


    def minimizeLCOEGreenfield(self, time_years = 1.):
        optMdot = self.full_system_solver.minimizeLCOEGreenfield(time_years = time_years)
        if optMdot:
            output = self.full_system.gatherOutput()
            output.optMdot = optMdot
            return output
        else:
            raise Exception('No minimum found!')


    def maximizePower(self, time_years = 1.):
        optMdot = self.full_system_solver.maximizePower(time_years = time_years)
        if optMdot:
            output = self.full_system.gatherOutput()
            output.optMdot = optMdot
            return output
        else:
            raise Exception('No maximum found!')


    def solve(self, m_dot, time_years):
        self.full_system.solve(m_dot = m_dot, time_years = time_years)
        output = self.full_system.gatherOutput()
        return output
