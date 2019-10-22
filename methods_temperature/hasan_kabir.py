import math
import scipy.optimize as optimize
import constants as const
import conversion_factors as conver


class HasanKabirMethod(object):

    def __init__(self, well):
        self.well = well
        self.sine = self.sine()
        self.coordinate = None
        self.pressure = None

    def sine(self):
        well = self.well
        pipe_casing = well.pipe_casing
        angle_horizontal = pipe_casing.angle_horizontal
        angle_horizontal *= conver.deg_to_rad
        sine = math.sin(angle_horizontal)
        return sine

    def target_function(self, temperature):
        well = self.well
        gas = well.gas
        pressure = self.pressure
        heat_capacity_specific = gas.heat_capacity_specific(pressure, temperature)
        heat_capacity_specific *= conver.kJ_to_J
        parameter_A = well.parameter_A(pressure, temperature)
        rock = well.rock
        temperature_formation_initial = rock.temperature_formation_initial
        gradient_temperature_geothermal = rock.gradient_temperature_geothermal
        coordinate = self.coordinate
        sine = self.sine
        term1 = temperature_formation_initial - gradient_temperature_geothermal * coordinate * sine
        term2 = parameter_A * (1 - math.exp(-coordinate / parameter_A))
        constant_gravity = const.CONSTANT_GRAVITY
        constant_gravity *= conver.bar_to_Pa
        term3 = gradient_temperature_geothermal * sine - constant_gravity * sine / heat_capacity_specific
        temperature_calculation = term1 + term2 * term3
        error_temperature = temperature - temperature_calculation
        return error_temperature

    def compute(self, coordinate, pressure):
        self.coordinate = coordinate
        self.pressure = pressure
        target_function = self.target_function
        temperature = optimize.root_scalar(target_function, method='bisect', bracket=[273.15, 1e3]).root
        return temperature
