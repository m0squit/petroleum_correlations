import ConversionFactors as conver


class Rock(object):
    """Класс описывает осредненную горную породу вокруг скважины.

    Attributes:
        temperature_formation_initial (float): Начальная темпераутра пласта, K.
        density (float): Плотность горной породы, kg/m3.
        gradient_temperature_geothermal (float): Геотермический градиент горной породы, K/m.
        thermal_conductivity (float): Теплопроводность горной породы, kJ/m/day/K.
        heat_capacity_specific(float): Удельная изобарная теплоемкость гороной породы, kJ/kg/K.

    """

    def __init__(self,
                 temperature_formation_initial):

        self.temperature_formation_initial = temperature_formation_initial
        self.density = 2504
        self.gradient_temperature_geothermal = 0.03 * (1 / 1)
        self.thermal_conductivity = 2.42 * (conver.J_to_kJ / 1 / conver.s_to_day / 1)
        self.heat_capacity_specific = 1256 * (conver.J_to_kJ / 1 / 1)