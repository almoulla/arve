import math
import numpy as np


class inverted_gaussian:
    def inverted_gaussian(self, x: np.ndarray, *params: np.float64) -> np.ndarray:
        """Inverted Gaussian.

        :param x: RV array
        :param params: tuple with continuum, contrast, RV and FWHM
        :type params: tuple of floats
        :return: inverted Gaussian evaluated at x
        """
        # unpack parameters
        continuum, contrast, RV, FWHM = params

        scale_factor = 1 / (2 * math.sqrt(2 * math.log(2)))  # incorrect log?
        # rename parameters
        C = continuum
        a = contrast
        b = RV
        # scale FWHM into sigma
        c = np.float64(FWHM * scale_factor)

        return C - np.exp(-(((x - b) / c) ** 2) / 2) * a
