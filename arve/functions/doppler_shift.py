import numpy as np


class doppler_shift:
    def doppler_shift(self, wave: np.ndarray, v: np.float64) -> np.ndarray:
        """Doppler shift.

        :param wave: wavelengths
        :param v: velocity in km/s
        :return: Doppler-shifted wavelengths
        """
        # vacuum speed of light
        c: float = 2.99792458e5  # [km/s]

        gamma_factor = (1 + v / c) / (1 - v / c)
        return wave * (gamma_factor) ** (1 / 2)
