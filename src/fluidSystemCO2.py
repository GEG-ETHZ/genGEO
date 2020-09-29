
from src.parasiticPowerFractionCoolingTower import parasiticPowerFractionCoolingTower
from src.powerPlantOutput import PowerPlantEnergyOutput

from utils.fluidStates import FluidState
from utils.fluidStateFromPT import FluidStateFromPT
from utils.fluidStateFromPh import FluidStateFromPh

class FluidSystemCO2Output(object):
    """FluidSystemCO2Output."""
    pass

class FluidSystemCO2(object):
    """FluidSystemCO2 provides methods to compute the water fluid cycle."""

    def __init__(self, T_ambient_C, dT_approach, eta_pump, eta_turbine, coolingMode):
        self.fluid = 'CO2'
        self.injection_well = None
        self.reservoir = None
        self.production_well = None

        self.T_ambient_C = T_ambient_C
        self.dT_approach = dT_approach
        self.eta_pump = eta_pump
        self.eta_turbine = eta_turbine
        self.coolingMode = coolingMode

    def solve(self, m_dot, time_years):

        self.pp_output = PowerPlantEnergyOutput()

        # Find condensation pressure
        T_condensation = self.T_ambient_C + self.dT_approach
        P_condensation = FluidState.getPFromTQ(T_condensation, 0, self.fluid) + 50e3
        dP_pump = 0

        dP_downhole_threshold = 1e4
        dP_downhole = dP_downhole_threshold
        while dP_downhole >= dP_downhole_threshold:

            # Find Injection Conditions
            P_pump_inlet = P_condensation
            P_pump_outlet = P_pump_inlet + dP_pump
            if P_pump_outlet < 7.38e6 and P_pump_outlet > 7.37e6:
                print('Total_analytic_system_co2:Manually adjusting pressure from ',
                        num2str(P_pump_outlet/1e6),
                        ' MPa to 7.37 MPa to avoid CoolProp CO2 critical point convergence issues.')
                P_pump_outlet = 7.37e6
            T_pump_inlet = T_condensation
            pump_inlet_state = FluidStateFromPT(P_pump_inlet, T_pump_inlet, self.fluid)
            h_pump_outletS = FluidState.getHFromPS(P_pump_outlet, pump_inlet_state.S_JK(), self.fluid)
            h_pump_outlet = pump_inlet_state.h_Jkg() + (h_pump_outletS - pump_inlet_state.h_Jkg()) / self.eta_pump

            self.pp_output.w_pump = -1 * (h_pump_outlet - pump_inlet_state.h_Jkg())

            surface_injection_state = FluidStateFromPh(P_pump_outlet, h_pump_outlet, self.fluid)

            injection_well_state    = self.injection_well.solve(surface_injection_state, m_dot, time_years)

            reservoir_state         = self.reservoir.solve(injection_well_state, m_dot, time_years)

            # find downhole pressure difference (negative means
            # overpressure
            dP_downhole = self.reservoir.P_reservoir - reservoir_state.P_Pa()
            #adjust pump if not within threshold
            if dP_downhole >= dP_downhole_threshold:
                dP_pump = dP_pump + 0.5 * dP_downhole

        if reservoir_state.P_Pa() >= self.reservoir.P_reservoir_max:
            raise ValueError('FluidSystemCO2:ExceedsMaxReservoirPressure - '
                        'Exceeds Max Reservoir Pressure of %.3f MPa!'%(self.reservoir.P_reservoir_max/1e6))

        production_well_state  = self.production_well.solve(reservoir_state, m_dot, time_years)

        # Calculate Turbine Power
        h_turbine_outS = FluidState.getHFromPS(P_condensation, production_well_state.S_JK(), self.fluid)
        h_turbine_out = production_well_state.h_Jkg() - self.eta_turbine * (production_well_state.h_Jkg() - h_turbine_outS)

        self.pp_output.w_turbine = production_well_state.h_Jkg() - h_turbine_out
        if self.pp_output.w_turbine < 0:
            raise ValueError('FluidSystemCO2:TurbinePowerNegative','Turbine Power is Negative')

        # heat rejection
        h_satVapor = FluidState.getHFromPQ(P_condensation, 1,  self.fluid)
        h_condensed = FluidState.getHFromPQ(P_condensation, 0,  self.fluid)
        if h_turbine_out > h_satVapor:
            # desuperheating needed
            self.pp_output.q_cooler = h_satVapor - h_turbine_out
            self.pp_output.q_condenser = h_condensed - h_satVapor
            T_turbine_out = FluidState.getTFromPh(P_condensation, h_turbine_out, self.fluid)
            dT_range = T_turbine_out - T_condensation
        else:
            # no desuperheating
            self.pp_output.q_cooler = 0
            self.pp_output.q_condenser = h_condensed - h_turbine_out
            dT_range = 0

        parasiticPowerFraction = parasiticPowerFractionCoolingTower(self.T_ambient_C, self.dT_approach, dT_range, self.coolingMode)
        self.pp_output.w_cooler = self.pp_output.q_cooler * parasiticPowerFraction('cooling')
        self.pp_output.w_condenser = self.pp_output.q_condenser * parasiticPowerFraction('condensing')

        self.pp_output.w_net = self.pp_output.w_turbine + self.pp_output.w_pump + self.pp_output.w_cooler + self.pp_output.w_condenser

        self.pp_output.dP_surface = production_well_state.P_Pa() - P_condensation

        self.pp_output.state_out = surface_injection_state
        return self.pp_output.state_out

    def gatherOutput(self):
        output = FluidSystemCO2Output()
        output.injection_well       = self.injection_well.gatherOutput()
        output.reservoir            = self.reservoir.gatherOutput()
        output.production_well      = self.production_well.gatherOutput()
        output.pp                   = self.pp_output
        return output
