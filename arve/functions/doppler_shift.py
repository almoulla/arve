class doppler_shift:

    def doppler_shift(self, wave:list, v:float) -> list:
        """Doppler shift.

        :param wave: wavelengths
        :type wave: list
        :param v: velocity in km/s
        :type v: float
        :return: Doppler-shifted wavelengths
        :rtype: list
        """

        # read constants
        c = self.constants["c"] # [km/s] speed of light in vacuum

        return wave*((1+v/c)/(1-v/c))**(1/2)