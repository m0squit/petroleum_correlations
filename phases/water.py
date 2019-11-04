import math
import conversion_factors as conver


class Water(object):
    """Класс описывает минерализованную пластовую воду.

    Attributes:
        salinity (float): Массовая доля соли в воде, unit fractions.

    """

    def __init__(self,
                 salinity):

        self.salinity = salinity

    @staticmethod
    def formation_volume_factor(pressure, temperature):
        p = pressure * conver.bar_to_psi
        t = temperature * 9 / 5 - 459.67
        term1 = -1.0001e-2 + 1.33391 * 1e-4 * t + 5.50654 * 1e-7 * (t ** 2)
        term2 = -1.95301e-9 * p * t - 1.72834e-13 * (p ** 2) * t - 3.58922e-7 * p - 2.25341e-10 * (p ** 2)
        formation_volume_factor = (1 + term1) * (1 + term2)
        return formation_volume_factor

    @classmethod
    def density_fresh(cls, pressure, temperature):
        p = pressure
        t = temperature
        term1 = -80 * t - 3.3 * (t ** 2) + 0.00175 * (t ** 3) + 489 * p
        term2 = -2 * t * p + 0.016 * (t ** 2) * p - 1.3e-5 * (t ** 3) * p - 0.333 * (p ** 2) - 0.002 * t * (p ** 2)
        density_fresh = 1 + 1e-6 * (term1 + term2)
        return density_fresh

    def density(self, pressure, temperature):
        p = pressure * conver.bar_to_MPa
        t = temperature - 273.15
        density_fresh = self.density_fresh(p, t)
        salinity = self.salinity
        term1 = 0.668 + 0.44 * salinity
        term2 = 1e-6 * (300 * p - 2400 * p * salinity + t * (80 + 3 * t - 3300 * salinity - 13 * p + 47 * p * salinity))
        density = density_fresh + salinity * (term1 + term2)
        density *= conver.g_per_cm3_to_kg_per_m3
        return density

    def viscosity_p0(self, temperature):
        t = temperature
        term1 = 20 - t
        term2 = 1.23780 * term1 - 1.303e-3 * (term1 ** 2) + 3.06e-6 * (term1 ** 3) + 2.550e-8 * (term1 ** 4)
        term3 = 96 + t
        term4 = term2 / term3
        viscosity_fresh_p0_t20 = 1.002
        viscosity_fresh_p0 = viscosity_fresh_p0_t20 * math.exp(term4)
        salinity = self.salinity
        salinity /= 0.05844
        a = 3.324e-2 * salinity + 3.624e-3 * (salinity ** 2) - 1.879e-4 * (salinity ** 3)
        b = -3.960e-2 * salinity + 1.020e-2 * (salinity ** 2) - 7.020e-4 * (salinity ** 3)
        term5 = a + b * term4
        viscosity_p0 = viscosity_fresh_p0 * math.exp(term5)
        return viscosity_p0

    def coefficient_viscosity(self, temperature):
        t = temperature
        term1 = -1.297 + 5.74e-2 * t - 6.97e-4 * (t ** 2) + 4.47e-6 * (t ** 3) - 1.05e-8 * (t ** 4)
        term2 = 0.545 + 2.8e-3 * t - term1
        term3 = 6.044 + 2.8e-3 * t + 3.6e-5 * (t ** 2)
        salinity = self.salinity
        salinity /= 0.05844
        term4 = salinity / term3
        term5 = 2.5 * term4 - 2 * (term4 ** 2) + 0.5 * (term4 ** 3)
        coefficient_viscosity = term2 * term5 + term1
        return coefficient_viscosity

    def viscosity(self, pressure, temperature):
        p = pressure * conver.bar_to_MPa
        t = temperature - 273.15
        viscosity_p0 = self.viscosity_p0(t)
        coefficient_viscosity = self.coefficient_viscosity(t)
        viscosity = viscosity_p0 * (1 + 0.001 * coefficient_viscosity * p)
        return viscosity
