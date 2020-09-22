import unittest

from src.semiAnalyticalWell import SemiAnalyticalWell
from src.porousReservoir import PorousReservoir
from src.subsurfaceComponents import DownHolePump
from src.oRCCycleTboil import ORCCycleTboil
from src.fluidSystemWater import FluidSystemWater
from src.fluidSystemWaterSolver import FluidSystemWaterSolver

from utils.globalProperties import GlobalSimulationProperties
from utils.globalConstants import globalConstants
from utils.fluidStateFromPT import FluidStateFromPT

from tests.testAssertion import testAssert


# define global methods to be used in this tests
gsp = GlobalSimulationProperties()
# alternative simulation properties for porous reservoir
gsp2 = GlobalSimulationProperties()
gsp2.k_rock = 2.1        #W/m/C
gsp2.rho_rock = 2300     #kg/m^3
gsp2.c_rock = 920.       #J/kg/C

inj_well = SemiAnalyticalWell(params = gsp,
                                    dz_total = -2500.,
                                    wellRadius = 0.205,
                                    wellMultiplier = 4.,
                                    fluid = 'Water',
                                    epsilon = 55 * 1e-6,
                                    dT_dz = 0.035,
                                    T_e_initial = 15.)

reservoir = PorousReservoir(params = gsp2,
                                well_spacing = 707.,
                                thickness = 100,
                                permeability = 1.0e-15 * 15000 / 100., # permeability = transmissivity / thickness
                                T_surface_rock = 15,
                                depth = 2500,
                                dT_dz = 0.035,
                                wellRadius = 0.205,
                                reservoirConfiguration = '5spot',
                                fluid = 'Water',
                                modelPressureTransient = False,
                                modelTemperatureDepletion = True)

prod_well1 = SemiAnalyticalWell(params = gsp,
                                    dz_total = 2000.,
                                    wellRadius = 0.205,
                                    wellMultiplier = 4.,
                                    fluid = 'Water',
                                    epsilon = 55 * 1e-6,
                                    dT_dz = 0.035,
                                    T_e_initial = 15. + 0.035 * 2500.)

prod_well2 = SemiAnalyticalWell(params = gsp,
                                    dz_total = 500.,
                                    wellRadius = 0.205,
                                    wellMultiplier = 4.,
                                    fluid = 'Water',
                                    epsilon = 55 * 1e-6,
                                    dT_dz = 0.035,
                                    T_e_initial = 15. + 0.035 * 500.)

pump = DownHolePump(well = prod_well2,
                    pump_depth = 500.,
                    max_pump_dP = 10.e6,
                    eta_pump = 0.75)

cycle = ORCCycleTboil(T_ambient_C = 15.,
                        dT_approach = 7.,
                        dT_pinch = 5.,
                        eta_pump = 0.9,
                        eta_turbine = 0.8,
                        coolingMode = 'Wet',
                        orcFluid = 'R245fa')

fluid_system = FluidSystemWater()
fluid_system.injection_well = inj_well
fluid_system.reservoir = reservoir
fluid_system.production_well1 = prod_well1
fluid_system.pump = pump
fluid_system.orc = cycle

initialState = FluidStateFromPT(1.e6, 60., fluid_system.fluid)

class FluidSystemWaterTest(unittest.TestCase):

    def testFluidSystemWaterMdot1(self):
        fluid_system_results = fluid_system.solve(initial_state = initialState,
                        m_dot = 1,
                        time_years = 1)

        self.assertTrue(*testAssert(fluid_system_results.P_Pa(), 1.352117133629285e+06, 'test_subsurface1_pressure'))
        self.assertTrue(*testAssert(fluid_system_results.T_C(), 54.038702389379260, 'test_subsurface1_temp'))

    def testFluidSystemWaterSolverMdot1(self):
        solver = FluidSystemWaterSolver(fluid_system)
        solver.solve(m_dot = 1, time_years = 1)

        output = solver.fluid_system.gatherOutput()

        self.assertTrue(*testAssert(output.production_well2.end_P_Pa(), 1.400500841642725e+06, 'test_subsurface_solver1_pressure'))
        self.assertTrue(*testAssert(output.production_well2.end_T_C(), 81.255435963675800, 'test_subsurface_solver1_temp'))
        self.assertTrue(*testAssert(output.orc.w_net, 5.3517e3, 'test_subsurface_solver1_w_net'))

    def testFluidSystemWaterMdot40(self):
        fluid_system_results = fluid_system.solve(initial_state = initialState,
                        m_dot = 40,
                        time_years = 1)

        self.assertTrue(*testAssert(fluid_system_results.P_Pa(), 8.144555604792739e+06, 'test_subsurface2_pressure'))
        self.assertTrue(*testAssert(fluid_system_results.T_C(), 58.834943596593840, 'test_subsurface2_temp'))

    def testFluidSystemWaterSolverMdot40(self):
        solver = FluidSystemWaterSolver(fluid_system)
        solver.solve(m_dot = 40, time_years = 1)

        output = solver.fluid_system.gatherOutput()

        print(*testAssert(output.production_well2.end_P_Pa(), 8.144555268219999e+06, 'test_subsurface_solver2_pressure'))
        self.assertTrue(*testAssert(output.production_well2.end_T_C(), 100.3126791060898, 'test_subsurface_solver2_temp'))
        print(*testAssert(output.orc.w_net * 40., 1.5416e5, 'test_subsurface_solver2_w_net'))
