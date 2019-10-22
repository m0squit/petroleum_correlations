import math
import scipy.optimize as optimize
import constants as const
import conversion_factors as conver
from methods_pressure import average
from methods_pressure import fanning
from methods_temperature import constant
from methods_temperature import hasan_kabir

class Well(object):
    """Класс описывает газовую скважину.

    Attributes:
        pipe_casing (pipe): Обсадная труба скважины, object.
        pipe_production (pipe): Насосно-компрессорная труба скважины, object.
        rock (rock): Осредненная горная порода вокруг скважины, object.
        gas (gas): Газ скважины, object.
        time_work (float): Время работы скважины с момента последнего запуска, day.
        rate_standard (float): Дебит скважины в стандартных термобарических условиях, sm3/day.
        pressure_wellhead (float): Давление на устье скважины, barsa.

    """

    def __init__(self,
                 pipe_casing,
                 pipe_production,
                 rock,
                 gas,
                 time_work,
                 rate_standard,
                 pressure_wellhead):

        self.pipe_casing = pipe_casing
        self.pipe_production = pipe_production
        self.rock = rock
        self.gas = gas
        self.time_work = time_work
        self.rate_standard = rate_standard
        self.pressure_wellhead = pressure_wellhead

    def number_reynolds(self, pressure, temperature):
        """Расчет числа Рейнольдса для потока газа в НКТ скважины.

        Note:
            Источник: https://en.wikipedia.org/wiki/Reynolds_number.

        Args:
            pressure (float): Давление газа в НКТ скважины, barsa.
            temperature (float): Температура газа в НКТ скважины, K.

        Returns:
            number_reynolds (float): Число Рейнольдса для потока газа в НКТ скважины
                в заданных термобарических условиях, dimensionless.

        """
        gas = self.gas
        density = gas.density(pressure, temperature)
        pipe_production = self.pipe_production
        diameter_inner = pipe_production.diameter_inner
        density_standard = gas.density(const.PRESSURE_STANDARD, const.TEMPERATURE_STANDARD)
        rate_standard = self.rate_standard
        rate = density_standard * rate_standard / density
        area_flow = math.pi * diameter_inner ** 2 / 4
        velocity = rate / area_flow
        velocity *= (1 / conver.day_to_s)
        viscosity = gas.viscosity(pressure, temperature)
        viscosity *= conver.cP_to_Pa_s
        number_reynolds = density * velocity * diameter_inner / viscosity
        return number_reynolds

    def number_nusselt(self, pressure, temperature):
        """Расчет числа Нуссельта для газа в НКТ скважины.

        Note:

        Args:

        Returns:
            number_nusselt (float): Число Нуссельта для газа в НКТ скважины
                в заданных термобарических условиях, dimensionless.

        """
        number_reynolds = self.number_reynolds(pressure, temperature)
        gas = self.gas
        number_prandtl = gas.number_prandtl(pressure, temperature)
        number_nusselt = 0.023 * number_reynolds ** 0.8 * number_prandtl ** 0.3
        return number_nusselt

    def number_grashof(self, pressure, temperature):
        """Расчет числа Грасгофа для газа в затрубном пространстве скважины.

        Note:
            Источник: https://en.wikipedia.org/wiki/Grashof_number.

        Args:

        Returns:
            number_grashof (float): Число Грасгофа для газа в затрубном пространстве сквжины
                в заданных термобарических условиях, dimensionless.

        """
        pipe_casing = self.pipe_casing
        radius_casing_inner = pipe_casing.diameter_inner / 2
        pipe_production = self.pipe_production
        radius_production_outer = pipe_production.diameter_outer / 2
        gas = self.gas
        density = gas.density(pressure, temperature)
        thermal_expansion = gas.thermal_expansion(pressure, temperature)
        viscosity = gas.viscosity(pressure, temperature)
        viscosity *= conver.cP_to_Pa_s
        term1 = (radius_casing_inner - radius_production_outer) ** 3 * density ** 2
        constant_gravity = const.CONSTANT_GRAVITY
        constant_gravity *= conver.bar_to_Pa
        term2 = thermal_expansion * constant_gravity * 3 / viscosity ** 2
        number_grashof = term1 * term2
        return number_grashof

    def coeff_heat_transfer_convective_forced(self, pressure, temperature):
        """Расчет коэффициента вынужденной конвективной теплопередачи для потока газа в НКТ.

        Note:
            Вынужденная (forced) конвективная теплопередача в потоке газа происходит,
                когда в процессе течения преобладают инерционные силы.
            Критерий: number_grashof / number_reynolds << 1.
            Источник: E.M.Al-Safran, J.P.Brill Applied Multiphase Flow in Pipes and Flow Assurance -
                Oil and Gas Production, стр. 145.

        Args:
            pressure (float): Давление газа, barsa.
            temperature (float): Температура газа, K.

        Returns:
            coeff_heat_transfer_convective_forced (float): Коэффициент вынужденной конвективной теплопередачи
                для потока газа в НКТ в заданных термобарических условиях, kJ/m2/day/K.

        """
        number_nusselt = self.number_nusselt(pressure, temperature)
        gas = self.gas
        thermal_conductivity = gas.thermal_conductivity(pressure, temperature)
        pipe_production = self.pipe_production
        diameter_inner = pipe_production.diameter_inner
        coeff_heat_transfer_convective_forced = number_nusselt * thermal_conductivity / diameter_inner
        return coeff_heat_transfer_convective_forced

    def coeff_heat_transfer_convective_natural(self, pressure, temperature):
        """Расчет коэффициента естественной конвективной теплопередачи для столба газа в затрубном пространстве.

        Note:
            Естественная (natural) конвективная теплопередача в столбе газа происходит в результате
                разной температуры между внешней стенкой НКТ и внутренней стенкой обсадной колонны.
                Более нагретый газ возле стенки НКТ становится легче и выталкивается наверх.
            Критерий: number_grashof / number_reynolds >> 1.
            Источник: E.M.Al-Safran, J.P.Brill Applied Multiphase Flow in Pipes and Flow Assurance -
                Oil and Gas Production, стр. 145.

        Args:
            pressure (float): Давление газа, barsa.
            temperature (float): Температура газа, K.

        Returns:
            coeff_heat_transfer_convective_natural (float): Коэффициент естественной конвективной теплопередачи
                для потока газа в затрубном пространстве в заданных термобарических условиях, kJ/m2/day/K.

        """
        number_grashof = self.number_grashof(pressure, temperature)
        gas = self.gas
        number_prandtl = gas.number_prandtl(pressure, temperature)
        thermal_conductivity = gas.thermal_conductivity(pressure, temperature)
        pipe_casing = self.pipe_casing
        radius_casing_inner = pipe_casing.diameter_inner / 2
        pipe_production = self.pipe_production
        radius_production_outer = pipe_production.diameter_outer / 2
        term1 = 0.049 * (number_grashof * number_prandtl) ** (1 / 3)
        term2 = number_prandtl ** 0.074 * thermal_conductivity
        term3 = radius_production_outer * math.log(radius_casing_inner / radius_production_outer)
        coeff_heat_transfer_convective_natural = term1 * term2 / term3
        coeff_heat_transfer_convective_natural *= 0.25
        return coeff_heat_transfer_convective_natural

    def time_transient(self):
        time_work = self.time_work
        rock = self.rock
        thermal_conductivity = rock.thermal_conductivity
        density = rock.density
        heat_capacity_specific = rock.heat_capacity_specific
        pipe_casing = self.pipe_casing
        cementing = pipe_casing.insulation
        radius_wellbore = cementing.diameter / 2
        time_transient = thermal_conductivity * time_work / (density * heat_capacity_specific * radius_wellbore ** 2)
        return time_transient

    @staticmethod
    def dissipation_heat_rock(time_transient):
        if time_transient <= 1.5:
            dissipation_heat_rock = 1.1281 * time_transient ** 0.5 * (1 - 0.3 * time_transient ** 0.5)
        else:
            dissipation_heat_rock = (0.4063 + 0.5 * math.log(time_transient)) * (1 + (0.6 / time_transient))
        return dissipation_heat_rock

    @staticmethod
    def resistance_thermal_convection(diameter_inner, coeff_heat_transfer_convective):
        radius_inner = diameter_inner / 2
        resistance_thermal_convection = 1 / (radius_inner * coeff_heat_transfer_convective)
        return resistance_thermal_convection

    def resistance_thermal(self, pressure, temperature):
        pipe_production = self.pipe_production
        diameter_inner_pipe_production = pipe_production.diameter_inner
        coeff_heat_transfer_convective_forced = self.coeff_heat_transfer_convective_forced(pressure, temperature)
        resistance_thermal = self.resistance_thermal_convection(diameter_inner_pipe_production, coeff_heat_transfer_convective_forced)
        resistance_thermal += pipe_production.resistance_thermal()

        pipe_casing = self.pipe_casing
        diameter_inner_pipe_casing = pipe_casing.diameter_inner
        coeff_heat_transfer_convective_natural = self.coeff_heat_transfer_convective_natural(pressure, temperature)
        resistance_thermal += self.resistance_thermal_convection(diameter_inner_pipe_casing, coeff_heat_transfer_convective_natural)
        resistance_thermal += pipe_casing.resistance_thermal()

        time_transient = self.time_transient()
        dissipation_heat_rock = self.dissipation_heat_rock(time_transient)
        rock = self.rock
        thermal_conductivity_rock = rock.thermal_conductivity
        resistance_thermal += dissipation_heat_rock / thermal_conductivity_rock
        return resistance_thermal

    def coeff_heat_transfer_overall(self, pressure, temperature):
        """Расчет суммарного коэффциента теплопередачи скважины.

        Returns:
            coeff_heat_transfer_overall (float): Cуммарный коэффциент теплопередачи
                в заданных термобарических условиях НКТ скважины, kJ/m2/day/K.

        """
        resistance_thermal = self.resistance_thermal(pressure, temperature)
        pipe_production = self.pipe_production
        radius_outer = pipe_production.diameter_outer / 2
        coeff_heat_transfer_overall = 1 / (radius_outer * resistance_thermal)
        return coeff_heat_transfer_overall

    def parameter_A(self, pressure, temperature):
        """Расчет параметра A.

        Note:
            Это обратная величина времени релаксации.

        Returns:
            parameter_A (float): Параметр A, m.

        """
        gas = self.gas
        heat_capacity_specific = gas.heat_capacity_specific(pressure, temperature)
        rate_standard = self.rate_standard
        density_standard = gas.density(const.PRESSURE_STANDARD, const.TEMPERATURE_STANDARD)
        rate_mass = density_standard * rate_standard
        coeff_heat_transfer_overall = self.coeff_heat_transfer_overall(pressure, temperature)
        pipe_production = self.pipe_production
        diameter_outer = pipe_production.diameter_outer
        parameter_A = heat_capacity_specific * rate_mass / (coeff_heat_transfer_overall * math.pi * diameter_outer)
        return parameter_A

    @staticmethod
    def select_method_pressure(name_method):
        methods_pressure = {'average': average.AverageMethod,
                            'fanning': fanning.FanningMethod}
        method_pressure = methods_pressure[name_method]
        return method_pressure

    @staticmethod
    def select_method_temperature(name_method):
        methods_temperature = {'constant': constant.ConstantMethod,
                               'hasan_kabir': hasan_kabir.HasanKabirMethod}
        method_temperature = methods_temperature[name_method]
        return method_temperature

    def compute_pressure_profile(self, method_pressure='average', method_temperature='constant', length_segment=10):
        """Расчет профиля давления по стволу скважины.

        Args:
            method_pressure (string): Название метода расчета давления.
            method_temperature (string): Название метода расчета температуры.
            length_segment (float): Длина сегмента колонны труб, m.

        """
        _method_pressure = self.select_method_pressure(method_pressure)
        _method_temperature = self.select_method_temperature(method_temperature)
        calculator_pressure = _method_pressure(self)
        calculator_temperature = _method_temperature(self)

        pipe_casing = self.pipe_casing
        length_pipe = pipe_casing.length
        coordinate_output = length_pipe
        pressure_output = self.pressure_wellhead
        temperature_output = calculator_temperature.compute(coordinate_output, pressure_output)

        profile = dict(coordinate=[], pressure=[], temperature=[])
        profile['coordinate'].append(coordinate_output)
        profile['pressure'].append(pressure_output)
        profile['temperature'].append(temperature_output)

        coordinate_input = coordinate_output - length_segment

        while coordinate_input != 0:

            def target_function(pressure_input):
                temperature_input = calculator_temperature.compute(coordinate_input, pressure_input)
                temperature_average = (temperature_output + temperature_input) / 2
                _pressure_input = calculator_pressure.compute(length_segment, pressure_output, temperature_average, )
                error_pressure_input = pressure_input - _pressure_input
                return error_pressure_input

            pressure_input = optimize.root_scalar(target_function, method='bisect', bracket=[pressure_output, 1e4]).root
            temperature_input = calculator_temperature.compute(coordinate_input, pressure_input)
            profile['coordinate'].append(coordinate_input)
            profile['pressure'].append(pressure_input)
            profile['temperature'].append(temperature_input)
            pressure_output = pressure_input
            temperature_output = temperature_input
            coordinate_input -= length_segment
        return profile
