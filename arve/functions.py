import math
from typing import Optional, TypeVar, Union

import numpy as np
import numpy.typing as npt

from arve import ARVE

TFunctions = TypeVar("TFunctions", bound="Functions")
RT = TypeVar("RT")


class Functions:
    """ARVE Functions base-class."""

    def __init__(self: TFunctions, arve: ARVE) -> None:
        self.arve = arve

    def convert_air_to_vac(
        self: TFunctions, wave_air: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
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

    def convert_vac_to_air(
        self: TFunctions, wave_vac: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
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

    def doppler_shift(
        self: TFunctions, wave: npt.NDArray[np.float64], v: np.float64
    ) -> npt.NDArray[np.float64]:
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
        self: TFunctions, x: npt.NDArray[np.float64], *params: np.float64
    ) -> npt.NDArray[np.float64]:
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

    def gls_periodogram(
        self: TFunctions,
        time: npt.NDArray[np.float64],
        val: npt.NDArray[np.float64],
        err: Optional[npt.NDArray[np.float64]] = None,
        ofac: int = 1,
        normalize: bool = True,
        win_func: bool = False,
    ) -> Union[
        tuple[
            npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]
        ],
        tuple[
            npt.NDArray[np.float64],
            npt.NDArray[np.float64],
            npt.NDArray[np.float64],
            npt.NDArray[np.float64],
            npt.NDArray[np.float64],
            float,
        ],
    ]:
        """Generalized Lomb-Scargle (GLS) periodogram.

        :param time: time values
        :param val: values
        :param err: errors, defaults to None
        :param ofac: over-factorization, defaults to 1
        :param normalize: normalization of periodogram, defaults to True
        :param win_func: return window function, defaults to False
        :return: frequency, power spectrum and phase of periodogram; if win_func is True, the frequency, power spectrum and phase of the window function are returned as well
        """
        # if not provided, set uncertainties to unity
        if err is None:
            err = np.ones(len(time))

        # frequencies, power spectrum and phases of data
        freq, ps, phi = _gls(time, val, err, ofac, normalize)

        # spectral window function
        if win_func:
            # central frequency
            freq_c = freq[len(freq) // 2]

            # sinusoid with unit amplitude at central frequency
            win_val = np.sin(2 * np.pi * freq_c * time)
            win_err = err

            # power spectrum of window function
            win_freq, win_ps, _ = _gls(time, win_val, win_err, ofac)

            # recenter window function frequencies
            win_freq -= freq_c

            # area of window function
            win_dfreq = np.mean(win_freq[1:] - win_freq[:-1])
            win_area = np.sum(win_ps) * win_dfreq

        # return periodogram parameters
        if win_func:
            return freq, ps, phi, win_freq, win_ps, win_area

        return freq, ps, phi


def _gls(
    time: npt.NDArray[np.float64],
    val: npt.NDArray[np.float64],
    err: npt.NDArray[np.float64],
    ofac: int,
    normalize: bool = True,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    # time span and steps
    Time = time[-1] - time[0]
    dtime = time[1:] - time[:-1]

    # linear and angular frequencies
    dfreq = 1 / (Time * ofac)
    freq = np.arange(1 / Time, 1 / (2 * np.median(dtime)), dfreq)
    omega = 2 * np.pi * freq

    # weights
    W = np.sum(1 / err**2)
    w = 1 / (W * err**2)

    # empty arrays for powers and phases
    N = len(freq)
    ps = np.empty(N)
    phi = np.empty(N)

    # loop through frequencies
    for i in range(N):
        # trigonometric terms
        arg = omega[i] * time
        cosarg = np.cos(arg)
        sinarg = np.sin(arg)

        # weighted sums
        Y = np.sum(w * val)
        C = np.sum(w * cosarg)
        S = np.sum(w * sinarg)

        # weighted sums of cross-terms
        YYhat = np.sum(w * val**2)
        YChat = np.sum(w * val * cosarg)
        YShat = np.sum(w * val * sinarg)
        CChat = np.sum(w * cosarg**2)
        SShat = np.sum(w * sinarg**2)
        CShat = np.sum(w * cosarg * sinarg)

        # differences of sums
        YY = YYhat - Y * Y
        YC = YChat - Y * C
        YS = YShat - Y * S
        CC = CChat - C * C
        SS = SShat - S * S
        CS = CShat - C * S

        # normalization
        D = CC * SS - CS**2

        # amplitudes
        a = (YC * SS - YS * CS) / D
        b = (YS * CC - YC * CS) / D

        # power spectrum
        if normalize:
            ps[i] = (SS * YC**2 + CC * YS**2 - 2 * CS * YC * YS) / (YY * D)
        if not normalize:
            ps[i] = a**2 + b**2

        # phases
        phi[i] = np.arctan2(a, b)

    # return frequencies, power spectrum and phases
    return freq, ps, phi
