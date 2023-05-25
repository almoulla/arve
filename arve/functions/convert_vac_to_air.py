class convert_vac_to_air:

    def convert_vac_to_air(self, wave_vac:list) -> list:
        """Convert wavelengths from vacuum to air.

        :param wave_vac: wavelengths in vacuum [Å]
        :type wave_vac: list
        :return: wavelengths in air [Å]
        :rtype: list
        """

        s = 1e4 / wave_vac
        n = 1 + 0.0000834254 + 0.02406147 / (130 - s**2) + 0.00015998 / (38.9 - s**2)
        
        wave_air = wave_vac / n

        return wave_air