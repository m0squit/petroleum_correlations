import math


class Pipe(object):
    """Класс описывает трубу любой секционной колонны труб скважины.

    Attributes:
        length (float): Длина трубы, m.
        diameter_outer (float): Наружный диаметр трубы, m.
        diameter_inner (float): Внутренний диаметр трубы, m.
        roughness_absolute (float): Абсолютная шероховатость стенок трубы, m.
        angle_horizontal (float): Угол отклонения трубы от горизонтали, deg.
        thermal_conductivity (float): Теплопроводность материала трубы, kJ/m/day/K.
        insulation (Insulation): Покрытие стенок трубы, object.

    """

    def __init__(self,
                 length,
                 diameter_outer,
                 diameter_inner,
                 roughness_absolute,
                 angle_horizontal,
                 thermal_conductivity,
                 insulation):

        self.length = length
        self.diameter_inner = diameter_inner
        self.diameter_outer = diameter_outer
        self.roughness_absolute = roughness_absolute
        self.angle_horizontal = angle_horizontal
        self.thermal_conductivity = thermal_conductivity
        self.insulation = insulation

    @staticmethod
    def resistance_thermal_conduction(diameter_outer, diameter_inner, thermal_conductivity):
        radius_outer = diameter_outer / 2
        radius_inner = diameter_inner / 2
        rt_conduction = math.log(radius_outer / radius_inner) / thermal_conductivity
        return rt_conduction

    def resistance_thermal(self):
        """Расчет термического сопротивления материала трубы с покрытием (если есть).

        Returns:
            resistance_thermal (float): Термическое сопротивление материала трубы с покрытием (если есть).

        """
        diameter_outer_pipe = self.diameter_outer
        diameter_inner_pipe = self.diameter_inner
        thermal_conductivity_pipe = self.thermal_conductivity
        resistance_thermal = self.resistance_thermal_conduction(diameter_outer_pipe, diameter_inner_pipe, thermal_conductivity_pipe)

        insulation = self.insulation
        if insulation is not None:
            diameter_outer_insulation = insulation.diameter
            thermal_conductivity_insulation = insulation.thermal_conductivity
            resistance_thermal += self.resistance_thermal_conduction(diameter_outer_insulation, diameter_outer_pipe, thermal_conductivity_insulation)
        return resistance_thermal

    def coeff_friction(self):
        diameter_inner = self.diameter_inner
        roughness_absolute = self.roughness_absolute
        coeff = 0.111 * (roughness_absolute * diameter_inner) ** 0.25
        return coeff