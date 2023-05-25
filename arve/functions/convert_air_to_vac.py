class convert_air_to_vac:

    def convert_air_to_vac(self, wave_air:list) -> list:
        """Convert wavelengths from air to vacuum.

        :param wave_air: wavelengths in air [Å]
        :type wave_air: list
        :return: wavelengths in vacuum [Å]
        :rtype: list
        """

        s = 1e4 / wave_air
        n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s**2) + 0.0001599740894897 / (38.92568793293 - s**2)
        
        wave_vac = wave_air * n

        return wave_vac