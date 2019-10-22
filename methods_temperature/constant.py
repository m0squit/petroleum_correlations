class ConstantMethod(object):

    def __init__(self, well):
        self.well = well

    @staticmethod
    def compute(*args):
        temperature_average = (349.9819 + 329.8312) / 2
        return temperature_average