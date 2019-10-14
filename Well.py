import math
import scipy.optimize as optimize
import Constants as const
import ConversionFactors as conver


class Well(object):
    """Класс описывает газовую скважину.

    Attributes:
        pipe_casing (Pipe): Обсадная труба скважины, object.
        pipe_production (Pipe): Насосно-компрессорная труба скважины, object.
        rock (Rock): Осредненная горная порода вокруг скважины, object.
        gas (Gas): Газ скважины, object.
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
            pressure (float): Давление газа, barsa.
            temperature (float): Температура газа, K.

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
        coeff_thermal_expansion = 0.004824
        coeff_thermal_expansion *= (1 / (1 + conver.degreeC_to_degreeK))
        viscosity = gas.viscosity(pressure, temperature)
        viscosity *= conver.cP_to_Pa_s
        term1 = (radius_casing_inner - radius_production_outer) ** 3 * density ** 2
        constant_gravity = const.CONSTANT_GRAVITY
        constant_gravity *= conver.bar_to_Pa
        term2 = coeff_thermal_expansion * constant_gravity * 3 / viscosity ** 2
        number_grashof = term1 * term2
        return number_grashof

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

    def temperature(self, coordinate, pressure):
        rock = self.rock
        temperature_formation_initial = rock.temperature_formation_initial
        gradient_temperature_geothermal = rock.gradient_temperature_geothermal
        pipe_casing = self.pipe_casing
        angle_horizontal = pipe_casing.angle_horizontal
        angle_horizontal *= conver.deg_to_rad
        sin = math.sin(angle_horizontal)

        def target_function(temperature):
            gas = self.gas
            heat_capacity_specific = gas.heat_capacity_specific(pressure, temperature)
            heat_capacity_specific *= conver.kJ_to_J
            parameter_A = self.parameter_A(pressure, temperature)
            term1 = temperature_formation_initial - gradient_temperature_geothermal * coordinate * sin
            term2 = parameter_A * (1 - math.exp(-coordinate / parameter_A))
            constant_gravity = const.CONSTANT_GRAVITY
            constant_gravity *= conver.bar_to_Pa
            term3 = gradient_temperature_geothermal * sin - constant_gravity * sin / heat_capacity_specific
            temperature_calculation = term1 + term2 * term3
            error_temperature = temperature - temperature_calculation
            return error_temperature

        temperature = optimize.root_scalar(target_function, method='bisect', bracket=[273.15, 1000]).root
        return temperature

    def simple_temperature(self, x, y):
        temperature_average = (349.9819 + 329.8312) / 2
        return temperature_average

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
        pipe_production = self.pipe_production
        roughness_absolute = pipe_production.roughness_absolute
        diameter_inner = pipe_production.diameter_inner
        roughness_relative = roughness_absolute / diameter_inner
        number_reynolds = self.number_reynolds(pressure, temperature)
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
        gas = self.gas
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
            pressure_difference_friction (float): Потеря давления флюида на трение в сегменте трубы, barsa.

        """
        friction_factor_fanning = self.friction_factor_fanning(pressure, temperature)
        gas = self.gas
        density = gas.density(pressure, temperature)
        density_standard = gas.density(const.PRESSURE_STANDARD, const.TEMPERATURE_STANDARD)
        rate_standard = self.rate_standard
        rate = density_standard * rate_standard / density
        pipe_production = self.pipe_production
        diameter_inner = pipe_production.diameter_inner
        area_flow = math.pi * diameter_inner ** 2 / 4
        velosity = rate / area_flow

        density *= (conver.kg_to_lb / conver.m_to_ft ** 3)
        velosity *= (conver.m_to_ft / conver.day_to_s)
        length_pipe_segment *= conver.m_to_ft
        constant_gravity = const.CONSTANT_GRAVITY
        constant_gravity *= (conver.bar_to_Pa * conver.m_to_ft)
        diameter_inner *= conver.m_to_ft

        term1 = 2 * friction_factor_fanning * density * velosity ** 2 * length_pipe_segment
        term2 = constant_gravity * diameter_inner
        pressure_loss_friction = term1 / term2
        pressure_loss_friction *= conver.psi_to_bar
        return pressure_loss_friction

    def fanning_correlation(self, pressure_output, temperature_average, length):

        def target_function(pressure_input):
            pressure_average = (pressure_output + pressure_input) / 2
            pressure_difference_hydrostatic = self.pressure_difference_hydrostatic(length, pressure_average, temperature_average)
            pressure_loss_friction = self.pressure_difference_friction(length, pressure_average, temperature_average)
            _pressure_input = pressure_output + (pressure_difference_hydrostatic + pressure_loss_friction)
            error_pressure_input = pressure_input - _pressure_input
            return error_pressure_input

        pressure_input = optimize.root_scalar(target_function, method='bisect', bracket=[pressure_output, 10000]).root
        return pressure_input

    def simple_correlation(self, pressure_output, temperature_average, length):
        gas = self.gas
        pipe_production = self.pipe_production

        density_relative = gas.density_relative
        angle_vertical = 90 - pipe_production.angle_horizontal
        angle_vertical *= conver.deg_to_rad
        cos = math.cos(angle_vertical)
        coeff_friction = pipe_production.coeff_friction()
        diameter_inner = pipe_production.diameter_inner

        def target_function(pressure_input):
            pressure_average = (pressure_output + pressure_input) / 2
            compressibility_factor = gas.compressibility_factor(pressure_average, temperature_average)
            parameter_s = 0.0375 * density_relative * length * cos / (compressibility_factor * temperature_average)
            density_standard = gas.density(const.PRESSURE_STANDARD, const.TEMPERATURE_STANDARD)
            rate_standard = self.rate_standard
            density = gas.density(pressure_average, temperature_average)
            rate = density_standard * rate_standard / density
            term1 = 6.67 * 1e-4 * rate ** 2 * coeff_friction * temperature_average ** 2 * compressibility_factor ** 2
            term2 = (math.exp(parameter_s) - 1) / diameter_inner ** 5 * cos
            _pressure_input = pressure_output ** 2 * math.exp(parameter_s) + term1 * term2
            error_pressure_input = pressure_input - _pressure_input
            return error_pressure_input

        pressure_input = optimize.root_scalar(target_function, method='bisect', bracket=[pressure_output, 10000]).root
        return pressure_input

    def compute_pressure_temperature_profile(self, method_pressure="fanning",
                                             method_temperature="simple", length_pipe_segment=10):

        methods_pressure = {"simple": self.simple_correlation,
                            "fanning": self.fanning_correlation}

        methods_temperature = {"simple": self.simple_temperature,
                               "hasan_kabir": self.temperature}

        pipe_casing = self.pipe_casing
        length_pipe = pipe_casing.length
        count_pipe_segment = int(length_pipe / length_pipe_segment)

        pt = dict(coordinate=[None] * count_pipe_segment,
                  pressure=[None] * count_pipe_segment,
                  temperature=[None] * count_pipe_segment)

        pressure_output = self.pressure_wellhead
        coordinate_output = length_pipe
        coordinate_input = length_pipe - length_pipe_segment
        for index_segment in range(count_pipe_segment):
            coordinate_average = (coordinate_output + coordinate_input) / 2
            pt['coordinate'][index_segment] = coordinate_average
            temperature_output = methods_temperature[method_temperature](coordinate_output, pressure_output)

            def target_function(pressure_input):
                temperature_input = methods_temperature[method_temperature](coordinate_input, pressure_input)
                temperature_average = (temperature_output + temperature_input) / 2
                pt['temperature'][index_segment] = temperature_average
                _pressure_input = methods_pressure[method_pressure](pressure_output, temperature_average, length_pipe_segment)
                error_pressure_input = pressure_input - _pressure_input
                return error_pressure_input

            pressure_input = optimize.root_scalar(target_function, method='bisect', bracket=[pressure_output, 10000]).root
            pressure_average = (pressure_output + pressure_input) / 2
            pt['pressure'][index_segment] = pressure_average
            pressure_output = pressure_input
            coordinate_output = coordinate_input
            coordinate_input -= length_pipe_segment
        return pt
