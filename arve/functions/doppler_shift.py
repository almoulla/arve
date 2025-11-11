import numpy as np

class doppler_shift:

    def doppler_shift(
        self,
        wave : float | np.ndarray,
        v    : float             ,
        ) -> float | np.ndarray:
        """Doppler shift.

        Parameters
        ----------
        wave : float | np.ndarray
            wavelength(s)
        v : float
            velocity in km/s

        Returns
        -------
        float | np.ndarray
            Doppler-shifted wavelength(s)
        """

        # read constants
        c = self.constants["c"] # [km/s] speed of light in vacuum

        return wave*((1+v/c)/(1-v/c))**(1/2)