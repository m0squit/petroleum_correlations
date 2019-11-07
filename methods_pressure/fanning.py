import math
import scipy.optimize as optimize
import constants as const
import conversion_factors as conver


class FanningMethod(object):

    def __init__(self, well):
        self.well = well
        self.length = None
        self.pressure_output = None
        self.temperature_average = None

    def friction_factor_fanning(self, pressure, temperature):
        """Расчет коэффициента трения для сегмента трубы по Fanning.

        Note:
            Источник: http://fekete.com/SAN/TheoryAndEquations/PiperTheoryEquations/c-te-pressure.htm.

        Args:
            pressure (float): Давление флюида, barsa.
            temperature (float): Температура флюида, K.

        Returns:
            friction_factor_fanning (float): Коэффициент трения для сегмента трубы с флюидом по Fanning.

        """
        well = self.well
        pipe_production = well.pipe_production
        roughness_absolute = pipe_production.roughness_absolute
        diameter_inner = pipe_production.diameter_inner
        roughness_relative = roughness_absolute / diameter_inner
        number_reynolds = well.number_reynolds(pressure, temperature)
        term1 = 0.2698 * roughness_relative
        term2 = 5.0452 / number_reynolds
        term3 = math.log(0.3539 * roughness_relative ** 1.1098 + 5.8506 / number_reynolds ** 0.8981)
        friction_factor_fanning = 1 / (-4 * math.log(term1 - term2 * term3)) ** 2
        return friction_factor_fanning

    def pressure_difference_hydrostatic(self, length_pipe_segment, pressure, temperature):
        """Расчет изменения гидростатического давления в сегменте трубы.

        Note:
            Источник: http://fekete.com/SAN/TheoryAndEquations/PiperTheoryEquations/c-te-pressure.htm.

        Args:
            length_pipe_segment (float): Длина сегмента трубы, m.
            pressure (float): Среднее давление флюида по длине сегмента, barsa.
            temperature (float): Средняя температура флюида по длине сегмента, K.

        Returns:
            pressure_loss_hydrostatic (float): Изменение гидростатического давления газа в сегменте трубы, barsa.

        """
        constant_gravity = const.CONSTANT_GRAVITY
        constant_gravity *= (conver.bar_to_Pa * conver.m_to_ft)
        well = self.well
        gas = well.fluid.phases[0]
        density = gas.density(pressure, temperature)
        density *= (conver.kg_to_lb / conver.m_to_ft ** 3)
        length_pipe_segment *= conver.m_to_ft
        pressure_loss_hydrostatic = constant_gravity / (144 * constant_gravity) * density * length_pipe_segment
        pressure_loss_hydrostatic *= conver.psi_to_bar
        return pressure_loss_hydrostatic

    def pressure_difference_friction(self, length_pipe_segment, pressure, temperature):
        """Расчет потери давления на трение в сегменте трубы.

        Note:
            Используется корреляция Fanning.
            Источник: http://fekete.com/SAN/TheoryAndEquations/PiperTheoryEquations/c-te-pressure.htm.

        Args:
            length_pipe_segment (float): Длина сегмента трубы, m.
            pressure (float): Среднее давление флюида по длине сегмента, barsa.
            temperature (float): Средняя температура флюида по длине сегмента, K.

        Returns:
            __calc_pressure_difference_friction (float): Потеря давления флюида на трение в сегменте трубы, barsa.

        """
        friction_factor_fanning = self.friction_factor_fanning(pressure, temperature)
        well = self.well
        gas = well.fluid.phases[0]
        density = gas.density(pressure, temperature)
        density_standard = gas.density(const.PRESSURE_STANDARD, const.TEMPERATURE_STANDARD)
        rate_standard = well.fluid.rates_standard[gas]
        rate = density_standard * rate_standard / density
        pipe_production = well.pipe_production
        diameter_inner = pipe_production.diameter_inner
        area_flow = math.pi * diameter_inner ** 2 / 4
        velocity = rate / area_flow

        density *= (conver.kg_to_lb / conver.m_to_ft ** 3)
        velocity *= (conver.m_to_ft / conver.day_to_s)
        length_pipe_segment *= conver.m_to_ft
        constant_gravity = const.CONSTANT_GRAVITY
        constant_gravity *= (conver.bar_to_Pa * conver.m_to_ft)
        diameter_inner *= conver.m_to_ft

        term1 = 2 * friction_factor_fanning * density * velocity ** 2 * length_pipe_segment
        term2 = constant_gravity * diameter_inner
        pressure_loss_friction = term1 / term2
        pressure_loss_friction *= conver.psi_to_bar
        return pressure_loss_friction

    def target_function(self, pressure_input):
        pressure_output = self.pressure_output
        pressure_average = (pressure_output + pressure_input) / 2
        length = self.length
        temperature_average = self.temperature_average
        pressure_difference_hydrostatic = self.pressure_difference_hydrostatic(length, pressure_average, temperature_average)
        pressure_loss_friction = self.pressure_difference_friction(length, pressure_average, temperature_average)
        _pressure_input = pressure_output + (pressure_difference_hydrostatic + pressure_loss_friction)
        error_pressure_input = pressure_input - _pressure_input
        return error_pressure_input

    def compute(self, length, pressure_output, temperature_average):
        self.length = length
        self.pressure_output = pressure_output
        self.temperature_average = temperature_average
        target_function = self.target_function
        pressure_input = optimize.root_scalar(target_function, method='brenth', bracket=[pressure_output, 1e3]).root
        return pressure_input
