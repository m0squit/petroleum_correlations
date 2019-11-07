import math
import scipy.interpolate as inter
import conversion_factors as conver


class Fluid(object):
    """Класс описывает флюид скважины.

    Attributes:
        pipe (Pipe): Труба течения флюида, object.
        pressure (float): Давление флюида, barsa.
        temperature (float): Температура флюида, K.
        phases (dict): Фазы в составе продукции, object.
        rates_standard (dict): Дебит скважины по фазам в стандартных условиях, sm3/day.

    """

    def __init__(self,
                 gas,
                 water,
                 rate_gas_standard,
                 rate_water_standard):

        self.pipe = None
        self.pressure = None
        self.temperature = None
        self.phases = [gas, water]
        self.rates_standard = {gas: rate_gas_standard, water: rate_water_standard}

    @staticmethod
    def calc_tension_gas_water(pressure, temperature):
        """Расчет поверхностного натяжения между газом и водой.

        Note:
            Источник: http://fekete.com/SAN/TheoryAndEquations/PiperTheoryEquations/c-te-pressure.htm

        Args:
            pressure (float): Давление флюида, barsa.
            temperature (float): Температура флюида, K.

        Returns:
            tension_gas_water (float): Поверхностное натяжение между газом и водой
                в заданных термобарических условиях, dynes/cm.

        """
        p = pressure * conver.bar_to_psi
        t = temperature * 9 / 5 - 459.67
        tension_gas_water = 0
        t1 = 74
        t2 = 280
        term1 = 75 - 1.108 * p ** 0.349
        term2 = 53 - 0.1048 * p ** 0.637
        if t < t1:
            tension_gas_water = term1
        if t > t2:
            tension_gas_water = term2
        if t1 < t < t2:
            f = inter.interp1d([t1, t2], [term1, term2])
            tension_gas_water = f(t)
        tension_gas_water *= conver.dynes_per_cm_to_lbf_per_s2
        return tension_gas_water

    def velocity_mixture(self, pressure, temperature):
        velocities_superficial = self.__calc_velocities_superficial(pressure, temperature)
        velocity_mixture = sum(velocities_superficial.values())
        return velocity_mixture

    def ratio_velocities(self, pressure, temperature):
        phases = self.phases
        gas = phases[0]
        water = phases[1]
        velocities_superficial = self.__calc_velocities_superficial(pressure, temperature)
        return velocities_superficial[water] / velocities_superficial[gas]

    def volume_fractions_input(self, pressure, temperature):
        phases = self.phases
        volume_fractions_input = dict.fromkeys(phases)
        velocities_superficial = self.__calc_velocities_superficial(pressure, temperature)
        velocity_mixture = self.velocity_mixture(pressure, temperature)
        for phase in phases:
            volume_fractions_input[phase] = velocities_superficial[phase] / velocity_mixture
        return volume_fractions_input

    def density_no_slip(self, pressure, temperature):
        phases = self.phases
        density_no_slip = 0
        volume_fractions_input = self.volume_fractions_input(pressure, temperature)
        for phase in phases:
            density = phase.density(pressure, temperature)
            density *= conver.kg_per_m3_to_lb_per_ft3
            density_no_slip += density * volume_fractions_input[phase]
        return density_no_slip

    def viscosity_no_slip(self, pressure, temperature):
        phases = self.phases
        viscosity_no_slip = 0
        volume_fractions_input = self.volume_fractions_input(pressure, temperature)
        for phase in phases:
            viscosity_no_slip += phase.viscosity(pressure, temperature) * volume_fractions_input[phase]
        return viscosity_no_slip

    def density(self, pressure, temperature, volume_fraction_water):
        phases = self.phases
        gas = phases[0]
        water = phases[1]
        density_gas = gas.density(pressure, temperature)
        density_water = water.density(pressure, temperature)
        density = density_gas * (1 - volume_fraction_water) + density_water * volume_fraction_water
        return density

    def __calc_velocities_superficial(self, pressure, temperature):
        pipe = self.pipe
        diameter_inner = pipe.diameter_inner
        area_flow = math.pi * diameter_inner ** 2 / 4
        phases = self.phases
        velocities_superficial = dict.fromkeys(phases)
        for phase in phases:
            rate_standard = self.rates_standard[phase]
            formation_volume_factor = phase.formation_volume_factor(pressure, temperature)
            rate = rate_standard * formation_volume_factor
            velocity_superficial = rate / area_flow
            velocity_superficial *= (1 / conver.day_to_s)
            velocities_superficial[phase] = velocity_superficial
        return velocities_superficial
