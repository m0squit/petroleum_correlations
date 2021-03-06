class Insulation(object):
    """Класс описывает любое покрытие стенок (как наружных, так и внутренних) трубы.

    Примеры:
    - цементажное покрытие обсадной колонны скважины
    - асфальтено-парафиновое покрытие внутренних стенок трубы
    - специальное покрытие внешних стенок наземных труб

    Attributes:
    diameter (float): Диаметр кольца покрытия, m.
        Диаметр отсчитывается от центральной оси трубы.
        Если диаметр кольца больше наружного диаметра трубы, то покрытие - наружное.
        Если диаметр кольца меньше внутреннего диаметра трубы, то покрытие - внутреннее.
    thermal_conductivity (float): Теплопроводность материала покрытия, kJ/m/day/K.

    """

    def __init__(self,
                 diameter,
                 thermal_conductivity):

        self.diameter = diameter
        self.thermal_conductivity = thermal_conductivity