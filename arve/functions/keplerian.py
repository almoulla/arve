import numpy as np

class keplerian:

    def keplerian(
        self,
        t       : np.ndarray  ,
        *params : tuple[float],
        ) -> np.ndarray:
        """Keplerian.

        Parameters
        ----------
        t : np.ndarray
            time values
        params : tuple[float]
            tuple with period, RV semi-amplitude, phase and RV offset

        Returns
        -------
        np.ndarray
            Keplerian evaluated at t
        """
        
        # unpack parameters
        P, K, p, C = params

        return K*np.sin(2*np.pi/P*t + p) + C