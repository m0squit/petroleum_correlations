import math
import scipy.optimize as optimize
import constants as const
import conversion_factors as conver


class GrayMethod(object):
    """Класс описывает метод Gray для расчета давления флюида в начале трубы.

    Attributes:


    """

    def __init__(self, well):
        self.pipe = well.pipe_production
        self.fluid = well.fluid
        self.length = None
        self.pressure_output = None
        self.temperature_average = None

    def compute(self, length, pressure_output, temperature_average):
        self.fluid.pipe = self.pipe
        self.length = length
        self.pressure_output = pressure_output
        self.temperature_average = temperature_average
        target_function = self.__target_function
        pressure_input = optimize.root_scalar(target_function, method='bisect', bracket=[pressure_output, 10000]).root
        return pressure_input

    def __target_function(self, pressure_input):
        pressure_output = self.pressure_output
        pressure_average = (pressure_output + pressure_input) / 2
        length = self.length
        temperature_average = self.temperature_average
        pressure_difference_hydrostatic = self.__calc_pressure_difference_hydrostatic(length, pressure_average, temperature_average)
        pressure_loss_friction = self.__calc_pressure_difference_friction(length, pressure_average, temperature_average)
        _pressure_input = pressure_output + (pressure_difference_hydrostatic + pressure_loss_friction)
        error_pressure_input = pressure_input - _pressure_input
        return error_pressure_input

    def __calc_pressure_difference_hydrostatic(self, length_pipe_segment, pressure, temperature):
        constant_gravity = const.CONSTANT_GRAVITY
        constant_gravity *= (conver.bar_to_Pa * conver.m_to_ft)
        fluid = self.fluid
        volume_fraction_water = self.__calc_volume_fraction_water(pressure, temperature)
        density = fluid.density(pressure, temperature, volume_fraction_water)
        density *= (conver.kg_to_lb / conver.m_to_ft ** 3)
        length_pipe_segment *= conver.m_to_ft
        pressure_loss_hydrostatic = constant_gravity / (144 * constant_gravity) * density * length_pipe_segment
        pressure_loss_hydrostatic *= conver.psi_to_bar
        return pressure_loss_hydrostatic

    def __calc_pressure_difference_friction(self, length_pipe_segment, pressure, temperature):
        friction_factor_fanning = self.__calc_friction_factor_fanning(pressure, temperature)
        fluid = self.fluid
        density_no_slip = fluid.density_no_slip(pressure, temperature)
        density_no_slip *= (conver.kg_to_lb / conver.m_to_ft ** 3)
        pipe = self.pipe
        diameter_inner = pipe.diameter_inner
        diameter_inner *= conver.m_to_ft
        velocity_mixture = fluid.velocity_mixture(pressure, temperature)
        velocity_mixture *= (conver.m_to_ft / conver.day_to_s)
        length_pipe_segment *= conver.m_to_ft
        constant_gravity = const.CONSTANT_GRAVITY
        constant_gravity *= (conver.bar_to_Pa * conver.m_to_ft)
        term1 = 2 * friction_factor_fanning * density_no_slip * velocity_mixture ** 2 * length_pipe_segment
        term2 = 144 * constant_gravity * diameter_inner
        pressure_loss_friction = term1 / term2
        pressure_loss_friction *= conver.psi_to_bar
        return pressure_loss_friction

    def __calc_number_1(self, pressure, temperature):
        fluid = self.fluid
        density_no_slip = fluid.density_no_slip(pressure, temperature)
        velocity_mixture = fluid.velocity_mixture(pressure, temperature)
        constant_gravity = const.CONSTANT_GRAVITY
        constant_gravity *= (conver.bar_to_Pa * conver.m_to_ft)
        surface_tension = fluid.surface_tension_gas_water(pressure, temperature)
        density_gas = fluid.phases[0].density(pressure, temperature)
        density_gas *= conver.kg_per_m3_to_lb_per_ft3
        density_water = fluid.phases[1].density(pressure, temperature)
        density_water *= conver.kg_per_m3_to_lb_per_ft3
        term1 = (density_no_slip ** 2) * (velocity_mixture ** 4)
        term2 = constant_gravity * surface_tension * (density_water - density_gas)
        number_1 = term1 / term2
        return number_1

    def __calc_number_2(self, pressure, temperature):
        fluid = self.fluid
        pipe = self.pipe
        constant_gravity = const.CONSTANT_GRAVITY
        constant_gravity *= (conver.bar_to_Pa * conver.m_to_ft)
        diameter_inner = pipe.diameter_inner
        diameter_inner *= conver.m_to_ft
        density_gas = fluid.phases[0].density(pressure, temperature)
        density_gas *= conver.kg_per_m3_to_lb_per_ft3
        density_water = fluid.phases[1].density(pressure, temperature)
        density_water *= conver.kg_per_m3_to_lb_per_ft3
        surface_tension = fluid.surface_tension_gas_water(pressure, temperature)
        number_2 = constant_gravity * (diameter_inner ** 2) * (density_water - density_gas) / surface_tension
        return number_2

    def __calc_number_3(self, pressure, temperature):
        fluid = self.fluid
        ratio_velocities = fluid.ratio_velocities(pressure, temperature)
        number_3 = 0.0814 * (1 - 0.0554 * math.log(1 + 730 * ratio_velocities / (ratio_velocities + 1)))
        return number_3

    def __calc_volume_fraction_water(self, pressure, temperature):
        fluid = self.fluid
        phases = fluid.phases
        water = phases[1]
        volume_fractions_input = fluid.volume_fractions_input(pressure, temperature)
        volume_fractions_input_water = volume_fractions_input[water]
        number_1 = self.__calc_number_1(pressure, temperature)
        number_2 = self.__calc_number_2(pressure, temperature)
        number_3 = self.__calc_number_3(pressure, temperature)
        term1 = -2.314 * (number_1 * (1 + 205 / number_2)) ** number_3
        volume_fraction_water = 1 - (1 - volume_fractions_input_water) * (1 - math.exp(term1))
        return volume_fraction_water

    def __calc_roughness_effective(self, pressure, temperature):
        fluid = self.fluid
        ratio_velocities = fluid.ratio_velocities(pressure, temperature)
        surface_tension_gas_water = fluid.surface_tension_gas_water(pressure, temperature)
        density_no_slip = fluid.density_no_slip(pressure, temperature)
        velocity_mixture = fluid.velocity_mixture(pressure, temperature)
        roughness_0 = 28.5 * surface_tension_gas_water / (density_no_slip * velocity_mixture ** 2)
        if ratio_velocities >= 0.007:
            roughness_effective = roughness_0
        else:
            pipe = self.pipe
            roughness_absolute = pipe.roughness_absolute
            roughness_absolute *= conver.m_to_ft
            roughness_effective = roughness_absolute + ratio_velocities * (roughness_0 - roughness_absolute) / 0.007
        return roughness_effective

    def __calc_friction_factor_fanning(self, pressure, temperature):
        pipe = self.pipe
        roughness_effective = self.__calc_roughness_effective(pressure, temperature)
        diameter_inner = pipe.diameter_inner
        diameter_inner *= conver.m_to_ft
        roughness_relative = roughness_effective / diameter_inner
        number_reynolds = 1e7
        term1 = 0.2698 * roughness_relative
        term2 = 5.0452 / number_reynolds
        term3 = math.log(0.3539 * roughness_relative ** 1.1098 + 5.8506 / number_reynolds ** 0.8981)
        friction_factor_fanning = 1 / (-4 * math.log(term1 - term2 * term3)) ** 2
        return friction_factor_fanning
