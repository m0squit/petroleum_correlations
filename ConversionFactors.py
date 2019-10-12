import math


kg_to_lb = 1 / 0.45359237  # https://en.wikipedia.org/wiki/Pound_(mass)
kg_to_g = 1e3
m_to_ft = 1 / 0.3048  # https://en.wikipedia.org/wiki/Foot_(unit)
ft_to_m = 0.3048  # https://en.wikipedia.org/wiki/Foot_(unit)
s_to_day = 1 / 86400
day_to_s = 86400
hr_to_day = 1 / 24
psia_to_barsa = 0.068947573  # https://en.wikipedia.org/wiki/Pounds_per_square_inch
barsa_to_Pa = 1e5
degreeK_to_degreeR = 9 / 5  # https://en.wikipedia.org/wiki/Conversion_of_units_of_temperature#Comparison_of_temperature_scales
degreeF_to_degreeK = (1 + 459.67) * 5 / 9  # https://en.wikipedia.org/wiki/Conversion_of_units_of_temperature#Comparison_of_temperature_scales
degreeC_to_degreeK = 273.15
degree_to_radian = math.pi / 180
m3_to_cm3 = 1e2 ** 3
kg_M_to_mol = 1e3  # Eclipse Technical Description 2014 (eng) стр. 985
Btu_to_kJ = 1.054615 / 1e3  # https://en.wikipedia.org/wiki/British_thermal_unit
cP_to_Pa_s = 1e-3
J_to_kJ = 1e-3
kJ_to_J = 1e3
