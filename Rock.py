class Rock(object):
    """Класс описывает осредненную горную породу вокруг скважины.

    Attributes:
        temperature_formation_initial (float): Начальная темпераутра пласта, K.

    """

    def __init__(self,
                 temperature_formation_initial):
        self.temperature_formation_initial = temperature_formation_initial
        self.density = 2504
        self.gradient_temperature_geothermal = 0.03
        self.thermal_conductivity = 2.42 * 0.001 * 86400 / (1 + 273.15)
        self.heat_capacity_specific = 1256 * 0.001 / (1 + 273.15)