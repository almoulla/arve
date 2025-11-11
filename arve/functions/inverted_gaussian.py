import numpy as np

class inverted_gaussian:

    def inverted_gaussian(
        self,
        x       : np.ndarray  ,
        *params : tuple[float],
        ) -> np.ndarray:
        """Inverted Gaussian.

        Parameters
        ----------
        x : np.ndarray
            abscissa values (e.g., RVs for CCF, or wavelengths for spectral lines)
        params : tuple[float]
            tuple with continuum, contrast, center and FWHM

        Returns
        -------
        np.ndarray
            inverted Gaussian evaluated at x
        """
        
        # unpack parameters
        continuum, contrast, center, fwhm = params
        
        # rename parameters
        C = continuum
        a = contrast
        b = center
        c = fwhm
        
        # scale FWHM into sigma
        c /= 2*np.sqrt(2*np.log(2))

        return C*(1-a*np.exp(-((x-b)/c)**2/2))