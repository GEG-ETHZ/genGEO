import unittest
import numpy as np

from src.semiAnalyticalWell import SemiAnalyticalWell
from src.porousReservoir import PorousReservoir
from src.fluidSystemWater import FluidSystemWater
from src.surfacePlantComponents import DownHolePump
from src.oRCCycleTboil import ORCCycleTboil

from utils.globalProperties import GlobalSimulationProperties
from utils.globalConstants import globalConstants
from utils.fluidStateFromPT import FluidStateFromPT

from tests.testAssertion import testAssert


# define global methods to be used in this tests
gpp = GlobalSimulationProperties()
gpp.k_rock = 2.1        #W/m/C
gpp.rho_rock = 2300     #kg/m^3
gpp.c_rock = 920.       #J/kg/C


inj_well = SemiAnalyticalWell(params = gpp,
                                    dz_total = -2500.,
                                    wellRadius = 0.205,
                                    fluid = 'Water',
                                    epsilon = 55 * 1e-6,
                                    dT_dz = 0.035,
                                    T_e_initial = 15.)

reservoir = PorousReservoir(
            params = gpp,
            well_spacing = 707.,
            thickness = 100,
            permeability = 1.5e-11,
            T_surface_rock = 15,
            depth = 2500,
            dT_dz = 0.035,
            wellRadius = 0.205,
            reservoirConfiguration = '5spot',
            fluid = 'Water',
            modelPressureTransient = False,
            modelTemperatureDepletion = False)

prod_well1 = SemiAnalyticalWell(params = gpp,
                                    dz_total = 2000.,
                                    wellRadius = 0.205,
                                    fluid = 'Water',
                                    epsilon = 55 * 1e-6,
                                    dT_dz = -0.035,
                                    T_e_initial = reservoir.reservoirT)

prod_well2 = SemiAnalyticalWell(params = gpp,
                                    dz_total = 500.,
                                    wellRadius = 0.205,
                                    fluid = 'Water',
                                    epsilon = 55 * 1e-6,
                                    dT_dz = -0.035,
                                    T_e_initial = reservoir.reservoirT)

pump = DownHolePump(well = prod_well2,
                    pump_depth = 500,
                    max_pump_dP = 10e6)

cycle = ORCCycleTboil(T_ambient_C = 15.,
                        dT_approach = 7.,
                        dT_pinch = 5.,
                        eta_pump = 0.9,
                        eta_turbine = 0.8,
                        coolingMode = 'Wet',
                        orcFluid = 'R245fa')

class FluidSystemWaterTest(unittest.TestCase):

    def testFluidSystemWater(self):

        fluid_system = FluidSystemWater()
        fluid_system.injection_well = inj_well
        fluid_system.reservoir = reservoir
        fluid_system.production_well1 = prod_well1
        fluid_system.pump = pump
        fluid_system.oRC = cycle

        initialState = FluidStateFromPT(1.e6, 15., fluid_system.fluid)
        results = fluid_system.solve(initial_state = initialState,
                        m_dot = 100,
                        time_years = 10)

        output = fluid_system.gatherOutput()

        # print(results)

        # print(output.injection_well.end_P_Pa())
        # print(output.injection_well.end_T_C())
        # print(output.injection_well.end_h_Jkg())
        #
        # fs = output.reservoir.finalState()
        # fp = output.reservoir.getPressure()
        # print(fs.S_JK, fp, output.reservoir.heat)
        #
        # print(output.production_well1.end_P_Pa())
        # print(output.production_well1.end_T_C())
        # print(output.production_well1.end_h_Jkg())
