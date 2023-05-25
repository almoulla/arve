import math
from typing import Callable, Type, TypeVar

import numpy as np

from arve import ARVE
from functions import gls_periodogram

TFunctions = TypeVar("TFunctions", bound="Functions")
RT = TypeVar("RT")


def add_methods(functions: list[Callable[..., RT]]) -> Callable[..., Type["Functions"]]:
    """Add methods to the base class."""

    def decorator(cls: Type[Functions]) -> Type[Functions]:
        for function in functions:
            setattr(cls, function.__name__, function)
        return cls

    return decorator


@add_methods(gls_periodogram)
class Functions:
    """ARVE Functions base-class."""

    def __init__(self: TFunctions, arve: ARVE) -> None:
        self.arve = arve

    def convert_air_to_vac(self: TFunctions, wave_air: np.ndarray) -> np.ndarray:
        """Convert wavelengths from air to vacuum.

        :param wave_air: wavelengths in air [Å]
        :return: wavelengths in vacuum [Å]
        """
        s = 1e4 / wave_air
        n = (
            1
            + 0.00008336624212083
            + 0.02408926869968 / (130.1065924522 - s**2)
            + 0.0001599740894897 / (38.92568793293 - s**2)
        )

        return wave_air * n

    def convert_vac_to_air(self: TFunctions, wave_vac: np.ndarray) -> np.ndarray:
        """Convert wavelengths from vacuum to air.

        :param wave_vac: wavelengths in vacuum [Å]
        :return: wavelengths in air [Å]
        """
        s = 1e4 / wave_vac
        n = (
            1
            + 0.0000834254
            + 0.02406147 / (130 - s**2)
            + 0.00015998 / (38.9 - s**2)
        )

        return wave_vac / n

    def doppler_shift(self: TFunctions, wave: np.ndarray, v: np.float64) -> np.ndarray:
        """Doppler shift.

        :param wave: wavelengths
        :param v: velocity in km/s
        :return: Doppler-shifted wavelengths
        """
        # vacuum speed of light
        c: float = 2.99792458e5  # [km/s]

        gamma_factor = (1 + v / c) / (1 - v / c)
        return wave * (gamma_factor) ** (1 / 2)

    def inverted_gaussian(
        self: TFunctions, x: np.ndarray, *params: np.float64
    ) -> np.ndarray:
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

        # how the f to type this correctly?
        return C - np.exp(-(((x - b) / c) ** 2) / 2) * a

    def sptype_to_num(self: TFunctions, sptype: str) -> int:
        """Spectral type to number.

        :param sptype: spectral type
        :type sptype: str
        :return: spectral type represented as a number
        :rtype: int
        """
        return "OBAFGKM".index(sptype[0]) * 10 + int(sptype[1])
