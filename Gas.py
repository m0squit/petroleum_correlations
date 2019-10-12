import math
import Constants as const
import ConversionFactors as conver


class Gas(object):
    """Класс описывает углеводородный газ.

    Attributes:
        pressure_pseudocritical (float): Псевдокритическое давление газа, barsa.
        temperature_pseudocritical (float): Псевдокритическая температура газа, K.
        density_relative (float): Относительная плотность газа по воздуху в стандартных условиях, dimensionless.
        molar_mass (float): Молярная масса газа, kg/kg-M.

    """

    def __init__(self,
                 pressure_pseudocritical,
                 temperature_pseudocritical,
                 density_relative):

        self.pressure_pseudocritical = pressure_pseudocritical
        self.temperature_pseudocritical = temperature_pseudocritical
        self.density_relative = density_relative
        self.molar_mass = const.MOLAR_MASS_AIR * density_relative

    def pressure_reduced(self, pressure):
        """Расчет приведенного давления газа.

        Args:
            pressure(float): Давление газа, barsa.

        Returns:
            pressure_reduced (float): Давление газа, приведенное к псевдокритическому давлению, dimensionless.

        """
        pressure_pseudocritical = self.pressure_pseudocritical
        pressure_reduced = pressure / pressure_pseudocritical
        return pressure_reduced

    def temperature_reduced(self, temperature):
        """Расчет приведенной температуры газа.

        Args:
            temperature (float): Температура газа, K.

        Returns:
            temperature_reduced (float): Температура газа, приведенная к псевдокритической температуре, dimensionless.

        """
        temperature_pseudocritical = self.temperature_pseudocritical
        temperature_reduced = temperature / temperature_pseudocritical
        return temperature_reduced

    def compressibility_factor(self, pressure, temperature):
        """Расчет коэффициента сверхсжимаемости газа.

        Note:
            Используется аналитическое выражение для коэффициента сверхсжимаемости сухого углеводородного газа.
            Источник: Методичка, стр. 10.

        Args:
            pressure (float): Давление газа, barsa.
            temperature (float): Темпераутра газа, K.

        Returns:
            compressibility_factor (float): Коэффициент сверхсжимемости газа
                в заданных термобарических условиях, unit fraction.

        """
        pressure_reduced = self.pressure_reduced(pressure)
        temperature_reduced = self.temperature_reduced(temperature)
        compressibility_factor = (0.4 * math.log10(temperature_reduced) + 0.73) ** pressure_reduced + 0.1 * pressure_reduced
        return compressibility_factor

    def density(self, pressure, temperature):
        """Расчет плотности газа.

        Note:
            Используестя обобщенное уравнение состояния для реального газа.

        Args:
            pressure (float): Давление газа, barsa.
            temperature (float): Темпераутра газа, K.

        Returns:
            density (float): Плотность газа в заданных термобарических условиях, kg/m3.

        """
        molar_mass = self.molar_mass
        compressibility_factor = self.compressibility_factor(pressure, temperature)
        density = pressure * molar_mass / (compressibility_factor * const.CONSTANT_GAS * temperature)
        return density

    def viscosity(self, pressure, temperature):
        """Расчет вязкости газа.

        Note:
            Используется корреляция Lee для природного газа.
            Источник: https://petrowiki.org/Gas_viscosity.

        Args:
            pressure (float): Давление газа, barsa.
            temperature (float): Темпераутра газа, K.

        Returns:
            viscosity (float): Вязкость газа в заданных термобарических условиях, cP.

        """
        density = self.density(pressure, temperature)
        density *= (conver.kg_to_g / conver.m3_to_cm3)
        molar_mass = self.molar_mass
        temperature *= conver.degreeK_to_degreeR
        k1 = (0.00094 + 2 * 1e-6 * molar_mass) * temperature ** 1.5 / (209 + 19 * molar_mass + temperature)
        x = 3.5 + 986 / temperature + 0.01 * molar_mass
        y = 2.4 - 0.2 * x
        viscosity = k1 * math.exp(x * density ** y)
        return viscosity

    def heat_capacity_specific(self, pressure, temperature):
        """Расчет удельной изобарной теплоемкости газа.

        Note:
            Используется корреляция L.A.Kareem, T.M.Iwalewa, J.E.Omeke для природного газа.
            Источник: L.A.Kareem, T.M.Iwalewa, J.E.Omeke Isobaric specific heat capacity of natural gas
                as a function of specific gravity, pressure and temperature стр. 4, 8, 9.

        Args:
            pressure (float): Давление газа, barsa.
            temperature (float): Темпераутра газа, K.

        Returns:
            heat_capacity_specific (float): Удельная изобарная теплоемкость газа
                в заданных термобарических условиях, kJ/kg/K.

        """
        density_relative = self.density_relative
        term1 = (-10.9602 * density_relative + 25.9033)
        term2 = (2.1517 * 1e-1 * density_relative - 6.8687 * 1e-2) * temperature
        term3 = (-1.3337 * 1e-4 * density_relative + 8.6387 * 1e-5) * temperature ** 2
        term4 = (3.1474 * 1e-8 * density_relative - 2.8396 * 1e-8) * temperature ** 3
        molar_heat_capacity_specific_ideal = term1 + term2 + term3 + term4

        a1 = 4.80828
        a2 = -4.01563
        a3 = -0.0700681
        a4 = 0.0567
        a5 = 2.36642
        a6 = -3.82421
        a7 = 7.71784
        r = 8.3145
        pressure_reduced = self.pressure_reduced(pressure)
        temperature_reduced = self.temperature_reduced(temperature)
        t = 1 / temperature_reduced
        term5 = (1 + (a1 * math.exp(a2 * (1 - t) ** 2) * pressure_reduced * t) ** 2)
        term6 = (a7 + a6 * (pressure_reduced * t) + a5 * (pressure_reduced * t) ** 2 + a4 * (pressure_reduced * t) ** 3)
        term7 = (a1 * math.exp(a2 * (1 - t) ** 2) * pressure_reduced * t) ** 2 * (a3 * (pressure_reduced * t) ** 6)
        term8 = term6 ** 3
        molar_heat_capacity_specific_residual = r * (term5 / term6 - term7 / term8)

        molar_heat_capacity_specific_real = molar_heat_capacity_specific_ideal + molar_heat_capacity_specific_residual
        molar_mass = self.molar_mass
        molar_mass *= (1 / conver.kg_M_to_mol)
        heat_capacity_specific = molar_heat_capacity_specific_real / molar_mass
        heat_capacity_specific /= 1e3
        return heat_capacity_specific

    def thermal_conductivity(self, pressure, temperature):
        """Расчет теплопроводности газа.

        Note:
            Используется корреляция A.Jarrahian, E.Heidaryan для природного газа.
            Источник: A.Jarrahian, E.Heidaryan A simple correlation to estimate natural gas thermal conductivity стр. 2.

        Args:
            pressure (float): Давление газа, barsa.
            temperature (float): Темпераутра газа, K.

        Returns:
            thermal_conductivity (float): Теплопроводность газа в заданных термобарических условиях, kJ/m/day/K.

        """
        pressure_reduced = self.pressure_reduced(pressure)
        temperature_reduced = self.temperature_reduced(temperature)

        a1 = 3.095251494612 * 1e-5
        a2 = -3.054731613002 * 1e-1
        a3 = 1.205296187262 * 1e-2
        a4 = -2.155542603544 * 1e-2
        density_relative = self.density_relative
        temperature *= conver.degreeK_to_degreeR
        term1 = (a1 * density_relative ** a2) * (temperature - 459.67)
        term2 = a3 + a4 * math.log(density_relative)
        k1_atm_uncorrected = term1 + term2

        mole_fraction_n2 = 0
        a5 = 1.695938319680 * 1e-2
        a6 = 1.983908703280 * 1e-3
        delta_k_n2 = mole_fraction_n2 * (a5 * math.log(density_relative) + a6)

        mole_fraction_co2 = 0
        a7 = 1.469572516483 * 1e-2
        a8 = -7.570807856000 * 1e-4
        delta_k_co2 = mole_fraction_co2 * (a7 * math.log(density_relative) + a8)

        k1_atm = k1_atm_uncorrected + delta_k_n2 + delta_k_co2

        a9 = 1.854452341597
        a10 = -1.275798197236 * 1e-3
        a11 = 1.925784814025 * 1e-1
        term3 = a9 / temperature_reduced ** 5
        term4 = pressure_reduced ** 4 / (temperature_reduced ** 20 + pressure_reduced ** 4)
        term5 = a10 * (pressure_reduced / temperature_reduced) ** 2
        term6 = a11 * (pressure_reduced / temperature_reduced)
        thermal_conductivity = k1_atm * (1 + term3 * term4 + term5 + term6)
        thermal_conductivity *= (conver.Btu_to_kJ/conver.ft_to_m/conver.hr_to_day) #/conver.degreeF_to_degreeK)
        return thermal_conductivity

    def number_prandtl(self, pressure, temperature):
        """Расчет числа Прандтля для газа.

        Note:
            Источник: https://en.wikipedia.org/wiki/Prandtl_number.

        Args:
            pressure (float): Давление газа, barsa.
            temperature (float): Температура газа, K.

        Returns:
            number_Prandtl (float): Число Прандтля для газа в заданных термобарических условиях, dimensionless.

        """
        heat_capacity_specific = self.heat_capacity_specific(pressure, temperature)
        viscosity = self.viscosity(pressure, temperature)
        viscosity *= conver.cP_to_Pa_s
        thermal_conductivity = self.thermal_conductivity(pressure, temperature)
        thermal_conductivity *= (1/1/conver.day_to_s/1)
        number_prandtl = heat_capacity_specific * viscosity / thermal_conductivity
        return number_prandtl

