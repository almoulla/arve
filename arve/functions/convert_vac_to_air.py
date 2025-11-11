import numpy as np

class convert_vac_to_air:

    def convert_vac_to_air(
        self,
        wave_vac : float | np.ndarray,
        ) -> float | np.ndarray:
        """Convert wavelengths from vacuum to air. Source: https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion

        Parameters
        ----------
        wave_vac : float | np.ndarray
            wavelength(s) in vacuum [Å]

        Returns
        -------
        float | np.ndarray
            wavelength(s) in air [Å]
        """

        s = 1e4 / wave_vac
        n = 1 + 0.0000834254 + 0.02406147 / (130 - s**2) + 0.00015998 / (38.9 - s**2)
        
        wave_air = wave_vac / n

        return wave_air