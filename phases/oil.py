class Oil(object):
    """Класс описывает углеводородную жидкость.

    """

    def __init__(self,
                 pressure_bubble,
                 density_relative):

        self.pressure_bubble = pressure_bubble
        self.density_relative = density_relative

    def density(self):
        pass

    def viscosity(self):
        pass
