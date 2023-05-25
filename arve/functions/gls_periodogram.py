import numpy as np
from typing import Optional, Union


class gls_periodogram:
    def gls_periodogram(
        self,
        time: np.ndarray,
        val: np.ndarray,
        err: Optional[np.ndarray] = None,
        ofac: int = 1,
        normalize: bool = True,
        win_func: bool = False,
    ) -> Union[
        tuple[np.ndarray, np.ndarray, np.ndarray],
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float],
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
    time: np.ndarray,
    val: np.ndarray,
    err: np.ndarray,
    ofac: int,
    normalize: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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
