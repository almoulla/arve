import numpy as np

class convert_air_to_vac:

    def convert_air_to_vac(
        self,
        wave_air : float | np.ndarray,
        ) -> float | np.ndarray:
        """Convert wavelengths from air to vacuum. Source: https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion

        Parameters
        ----------
        wave_air : float | np.ndarray
            wavelength(s) in air [Å]

        Returns
        -------
        float | np.ndarray
            wavelength(s) in vacuum [Å]
        """

        s = 1e4 / wave_air
        n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s**2) + 0.0001599740894897 / (38.92568793293 - s**2)
        
        wave_vac = wave_air * n

        return wave_vac