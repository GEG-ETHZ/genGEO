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

from src.coolingCondensingTower import CoolingCondensingTower
from src.powerPlantOutput import PowerPlantEnergyOutput

from utils.fluidState import FluidState
from utils.solver import Solver

from utils.frictionFactor import frictionFactor

class FluidSystemCO2Output(object):
    """FluidSystemCO2Output."""
    pass

class FluidSystemCO2(object):
    """FluidSystemCO2 provides methods to compute the water fluid cycle."""

    def __init__(self, params):

        self.params = params

        self.injection_well = None
        self.reservoir = None
        self.production_well = None

    def solve(self):

        results = FluidSystemCO2Output()
        results.pp = PowerPlantEnergyOutput()

        # Find condensation pressure
        T_condensation = self.params.T_ambient_C + self.params.dT_approach
        P_condensation = FluidState.getStateFromTQ(T_condensation, 0, self.params.working_fluid).P_Pa + 50e3
        dP_pump = 0

        dP_downhole_threshold = 1e3
        dP_downhole = np.nan
        dP_loops = 1
        dP_solver = Solver()
        while np.isnan(dP_downhole) or np.abs(dP_downhole) >= dP_downhole_threshold:

            # Find Injection Conditions
            P_pump_inlet = P_condensation

            # Only add pump differential if positive. If negative, add as a throttle at the bottom of the injection well
            if (dP_pump > 0):
                P_pump_outlet = P_pump_inlet + dP_pump
            else:
                P_pump_outlet = P_pump_inlet

            T_pump_inlet = T_condensation
            pump_inlet_state = FluidState.getStateFromPT(P_pump_inlet, T_pump_inlet, self.params.working_fluid)

            if dP_pump > 0:
                h_pump_outletS = FluidState.getStateFromPS(P_pump_outlet, pump_inlet_state.s_JK, self.params.working_fluid).h_Jkg
                h_pump_outlet = pump_inlet_state.h_Jkg + (h_pump_outletS - pump_inlet_state.h_Jkg) / self.params.eta_pump_co2
            else:
                h_pump_outlet = pump_inlet_state.h_Jkg

            results.pp.w_pump = -1 * (h_pump_outlet - pump_inlet_state.h_Jkg)

            surface_injection_state = FluidState.getStateFromPh(P_pump_outlet, h_pump_outlet, self.params.working_fluid)

            results.injection_well    = self.injection_well.solve(surface_injection_state)

            # if dP_pump is negative, this is a throttle after the injection well
            if (dP_pump < 0):
                results.injection_well_downhole_throttle = FluidState.getStateFromPh(
                    results.injection_well.state.P_Pa + dP_pump, 
                    results.injection_well.state.h_Jkg, 
                    self.params.working_fluid)
            else:
                results.injection_well_downhole_throttle = FluidState.getStateFromPh(
                    results.injection_well.state.P_Pa, 
                    results.injection_well.state.h_Jkg, 
                    self.params.working_fluid)

            results.reservoir = self.reservoir.solve(results.injection_well_downhole_throttle)

            # find downhole pressure difference (negative means
            # overpressure
            dP_downhole = self.params.P_reservoir() - results.reservoir.state.P_Pa
            dP_pump = dP_solver.addDataAndEstimate(dP_pump, dP_downhole)
            if np.isnan(dP_pump):
                dP_pump = 0.5 * dP_downhole

            if dP_loops > 10:
                print('GenGeo::Warning::FluidSystemCO2:dP_loops is large: %s'%dP_loops)
            dP_loops += 1

        if results.reservoir.state.P_Pa >= self.params.P_reservoir_max():
            raise Exception('GenGeo::FluidSystemCO2:ExceedsMaxReservoirPressure - '
                        'Exceeds Max Reservoir Pressure of %.3f MPa!'%(self.params.P_reservoir_max()/1e6))

        results.production_well  = self.production_well.solve(results.reservoir.state)

        # Subtract surface frictional losses between production wellhead and surface plant
        ff = frictionFactor(self.params.well_radius, results.production_well.state.P_Pa, results.production_well.state.h_Jkg,
            self.params.m_dot_IP, self.params.working_fluid, self.params.epsilon)
        if self.params.has_surface_gathering_system == True:
            dP_surfacePipes = ff * self.params.well_spacing / (self.params.well_radius*2)**5 * 8 * self.params.m_dot_IP**2 / results.production_well.state.rho_kgm3 / 3.14159**2
        else:
            dP_surfacePipes = 0
        
        results.surface_plant_inlet = FluidState.getStateFromPh(
            results.production_well.state.P_Pa - dP_surfacePipes,
            results.production_well.state.h_Jkg,
            self.params.working_fluid)

        # Calculate Turbine Power
        h_turbine_outS = FluidState.getStateFromPS(P_condensation, results.surface_plant_inlet.s_JK, self.params.working_fluid).h_Jkg
        h_turbine_out = results.surface_plant_inlet.h_Jkg - self.params.eta_turbine_co2 * (results.surface_plant_inlet.h_Jkg - h_turbine_outS)

        results.pp.w_turbine = results.surface_plant_inlet.h_Jkg - h_turbine_out
        if results.pp.w_turbine < 0:
            raise Exception('GenGeo::FluidSystemCO2:TurbinePowerNegative - Turbine Power is Negative')

        # heat rejection
        h_satVapor = FluidState.getStateFromPQ(P_condensation, 1,  self.params.working_fluid).h_Jkg
        h_condensed = FluidState.getStateFromPQ(P_condensation, 0,  self.params.working_fluid).h_Jkg
        if h_turbine_out > h_satVapor:
            # desuperheating needed
            results.pp.q_cooler = h_satVapor - h_turbine_out
            results.pp.q_condenser = h_condensed - h_satVapor
            T_turbine_out = FluidState.getStateFromPh(P_condensation, h_turbine_out, self.params.working_fluid).T_C
            dT_range = T_turbine_out - T_condensation
        else:
            # no desuperheating
            results.pp.q_cooler = 0
            results.pp.q_condenser = h_condensed - h_turbine_out
            dT_range = 0

        parasiticPowerFraction = CoolingCondensingTower.parasiticPowerFraction(self.params.T_ambient_C, self.params.dT_approach, dT_range, self.params.cooling_mode)
        results.pp.w_cooler = results.pp.q_cooler * parasiticPowerFraction('cooling')
        results.pp.w_condenser = results.pp.q_condenser * parasiticPowerFraction('condensing')

        results.pp.w_net = results.pp.w_turbine + results.pp.w_pump + results.pp.w_cooler + results.pp.w_condenser

        results.pp.dP_surface = results.surface_plant_inlet.P_Pa - P_condensation
        results.pp.dP_pump = dP_pump

        results.pp.state_out = surface_injection_state

        return results
