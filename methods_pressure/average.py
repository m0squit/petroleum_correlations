import math
import scipy.optimize as optimize
import constants as const
import conversion_factors as conver


class AverageMethod(object):
    """Класс описывает метод средних температуры и z-фактора газа
        для расчета давления газа в начале трубы.

    Attributes:
        length (float): Длина трубы, m.
        pressure_output (float): Давление газа в конце трубы, barsa.
        temperature_average (float): Средняя температура газа по длине трубы, K.
        well (Well): Скважина, в которой установлена труба, object.

    """

    def __init__(self, well):
        self.well = well
        self.cosine = self.cosine()
        self.length = None
        self.pressure_output = None
        self.temperature_average = None

    def cosine(self):
        well = self.well
        pipe_casing = well.pipe_casing
        angle_horizontal = pipe_casing.angle_horizontal
        angle_vertical = 90 - angle_horizontal
        angle_vertical *= conver.deg_to_rad
        cosine = math.cos(angle_vertical)
        return cosine

    def parameter_s(self, compressibility_factor):
        well = self.well
        gas = well.gas
        density_relative = gas.density_relative
        length = self.length
        length *= conver.m_to_ft
        cosine = self.cosine
        temperature_average = self.temperature_average
        temperature_average *= conver.degreeK_to_degreeR
        parameter_s = 0.0375 * density_relative * length * cosine / (compressibility_factor * temperature_average)
        return parameter_s

    def target_function(self, pressure_input):
        pressure_output = self.pressure_output
        pressure_output *= conver.bar_to_psi
        well = self.well
        gas = well.gas
        pressure_input *= conver.bar_to_psi
        pressure_average = (pressure_output + pressure_input) / 2
        temperature_average = self.temperature_average
        compressibility_factor = gas.compressibility_factor(pressure_average, temperature_average)
        parameter_s = self.parameter_s(compressibility_factor)
        pipe_production = well.pipe_production
        coeff_friction = pipe_production.coeff_friction()
        density_standard = gas.density(const.PRESSURE_STANDARD, const.TEMPERATURE_STANDARD)
        rate_standard = well.rate_standard
        density = gas.density(pressure_average, temperature_average)
        rate = density_standard * rate_standard / density
        rate *= (conver.m_to_ft ** 3 / 1e6)
        diameter_inner = pipe_production.diameter_inner
        diameter_inner *= conver.m_to_ft
        cosine = self.cosine
        term1 = 6.67 * 1e-4 * rate ** 2 * coeff_friction * temperature_average ** 2 * compressibility_factor ** 2
        term2 = (math.exp(parameter_s) - 1) / diameter_inner ** 5 * cosine
        _pressure_input = math.sqrt(pressure_output ** 2 * math.exp(parameter_s) + term1 * term2)
        _pressure_input *= conver.psi_to_bar
        pressure_input *= conver.psi_to_bar
        error_pressure_input = pressure_input - _pressure_input
        return error_pressure_input

    def compute(self, length, pressure_output, temperature_average):
        self.length = length
        self.pressure_output = pressure_output
        self.temperature_average = temperature_average
        target_function = self.target_function
        pressure_input = optimize.root_scalar(target_function, method='bisect', bracket=[pressure_output, 1e4]).root
        return pressure_input
