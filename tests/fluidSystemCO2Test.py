import unittest
import numpy as np

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

fluid_system = FluidSystemCO2(T_ambient_C = 15.,
                        dT_approach = 7.,
                        eta_pump = 0.9,
                        eta_turbine = 0.78,
                        coolingMode = 'Wet')

fluid_system.injection_well = SemiAnalyticalWell(params = gsp,
                                    dz_total = -2500.,
                                    wellRadius = 0.205,
                                    wellMultiplier = 4.,
                                    fluid = 'CO2',
                                    epsilon = 55 * 1e-6,
                                    dT_dz = 0.035,
                                    T_e_initial = 15.)

fluid_system.reservoir = PorousReservoir(params = gsp2,
                                well_spacing = 707.,
                                thickness = 100,
                                permeability = 1.0e-15 * 15000 / 100., # permeability = transmissivity / thickness
                                T_surface_rock = 15,
                                depth = 2500,
                                dT_dz = 0.035,
                                wellRadius = 0.205,
                                reservoirConfiguration = '5spot',
                                fluid = 'CO2',
                                modelPressureTransient = False,
                                modelTemperatureDepletion = True)

fluid_system.production_well = SemiAnalyticalWell(params = gsp,
                                    dz_total = 2500.,
                                    wellRadius = 0.205,
                                    wellMultiplier = 4.,
                                    fluid = 'CO2',
                                    epsilon = 55 * 1e-6,
                                    dT_dz = 0.035,
                                    T_e_initial = 15. + 0.035 * 2500.)


capital_cost_system = CapitalCostSystem(2019)
capital_cost_system.CapitalCost_SurfacePlant = CapitalCostSurfacePlantCPG(2019)
capital_cost_system.CapitalCost_SurfacePipe = CapitalCostSurfacePipes(N = 0)
capital_cost_system.CapitalCost_Production_Well = CapitalCostWell.cO2Baseline(well_length = 2500,
                                                                    well_diameter = 0.205 * 2,
                                                                    success_rate = 0.95,
                                                                    cost_year = 2019)
capital_cost_system.CapitalCost_Injection_Well = CapitalCostWell.cO2Baseline(well_length = 2500,
                                                                    well_diameter = 0.205 * 2,
                                                                    success_rate = 0.95,
                                                                    cost_year = 2019)
capital_cost_system.CapitalCost_Wellfield = CapitalCostWellField.cO2MonitoringBaseline(N=1,
                                                                    monitoring_well_length=2500,
                                                                    monitoring_well_diameter=0.216,
                                                                    cost_year=2019)
capital_cost_system.CapitalCost_Exploration = CapitalCostExploration.cO2Baseline(N=1,
                                                                    well_length = 2500,
                                                                    well_diameter = 0.205 * 2,
                                                                    success_rate = 0.95,
                                                                    cost_year = 2019)
capital_cost_system.lcoe_model = LCOESimple(F_OM = 0.045,
                                            discountRate = 0.096,
                                            Lifetime = 25,
                                            CapacityFactor = 0.9)

full_system = FullSystemCPG(fluid_system, capital_cost_system)

class FluidSystemCO2Test(unittest.TestCase):

    def testFluidSystemCO2Mdot10(self):
        full_system.solve(m_dot = 10, time_years = 1)

        output = full_system.gatherOutput()

        self.assertTrue(*testAssert(output.fluid_system_solver.pp.dP_surface, 6.3943e6, 'test_dP_surface'))
        self.assertTrue(*testAssert(output.fluid_system_solver.production_well.end_T_C(), 60.4516, 'test_T_prod_surface_C'))
        self.assertTrue(*testAssert(output.energy_results.W_net, 1.0283e5, 'test_W_net'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_brownfield, 1.5452e7, 'test_C_brownfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_greenfield, 3.8082e7, 'test_C_greenfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 7.2330e-4, 'test_LCOE_brownfield'))


    def testFluidSystemCO2Mdot80(self):
        full_system.solve(m_dot = 80, time_years = 1)

        output = full_system.gatherOutput()

        self.assertTrue(*testAssert(output.fluid_system_solver.pp.dP_surface, 5.65949586e6, 'test_dP_surface'))
        self.assertTrue(*testAssert(output.fluid_system_solver.production_well.end_T_C(), 59.0826, 'test_T_prod_surface_C'))
        self.assertTrue(*testAssert(output.energy_results.W_net, 5.02286e+05, 'test_W_net'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_brownfield, 2.98408178e+07, 'test_C_brownfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_greenfield, 5.2470e7, 'test_C_greenfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 2.8596677e-04, 'test_LCOE_brownfield'))

    def testFluidSystemCO2Mdot200(self):
        full_system.solve(m_dot = 200, time_years = 1)

        output = full_system.gatherOutput()

        self.assertTrue(*testAssert(output.fluid_system_solver.pp.dP_surface, 3.53347e+06, 'test_dP_surface'))
        self.assertTrue(*testAssert(output.fluid_system_solver.production_well.end_T_C(), 47.603999, 'test_T_prod_surface_C'))
        self.assertTrue(*testAssert(output.energy_results.W_net, -1.3588e6, 'test_W_net'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_brownfield, 5.2721e7, 'test_C_brownfield_N'))
        self.assertTrue(*testAssert(output.capital_cost_model.C_greenfield, 7.5351e7, 'test_C_greenfield_N'))
        self.assertTrue(np.isnan(output.capital_cost_model.LCOE_brownfield.LCOE), 'test_LCOE_brownfield')

    def testFluidSystemCO2SolverOptMdot(self):

        full_system_solver = FullSystemSolverMinLCOEBrownfield(full_system)

        optMdot = full_system_solver.solve(time_years = 1)

        output = full_system.gatherOutput()
        print(*testAssert(optMdot, 55.56, 'test_optMdot_solver_optMdot'))
        print(*testAssert(output.energy_results.W_net, 4.5795e5, 'test_optMdot_solver_w_net'))
        print(*testAssert(output.capital_cost_model.LCOE_brownfield.LCOE, 2.6158e-4, 'test_optMdot_solver_LCOE_brownfield'))
