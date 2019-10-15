import Rock as r
import Gas as g
import Insulation as i
import Pipe as p
import Well as w


rock = r.Rock(temperature_formation_initial=350)
gas = g.Gas(pressure_pseudocritical=46, temperature_pseudocritical=190.5, density_relative=0.56)
cementing = i.Insulation(diameter=0.25, thermal_conductivity=0.7)

pipe_casing = p.Pipe(length=1000,
                     diameter_outer=0.2,
                     diameter_inner=0.15,
                     roughness_absolute=None,
                     angle_horizontal=90,
                     thermal_conductivity=0.7,
                     insulation=cementing)

pipe_production = p.Pipe(length=1000,
                         diameter_outer=0.125,
                         diameter_inner=0.10,
                         roughness_absolute=1e-5,
                         angle_horizontal=90,
                         thermal_conductivity=0.7,
                         insulation=None)

well = w.Well(pipe_casing=pipe_casing,
              pipe_production=pipe_production,
              rock=rock,
              gas=gas,
              time_work=31,
              rate_standard=1000,
              pressure_wellhead=100)

a = well.compute_pressure_temperature_profile(method_temperature="hasan_kabir")
print(f"a = {a}")
