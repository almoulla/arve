import numpy as np

class gls_periodogram:

    def gls_periodogram(self, time:list, val:list, err:list=None, ofac:int=1, normalize:bool=True, win_func:bool=False) -> tuple:
        """Generalized Lomb-Scargle (GLS) periodogram.

        :param time: time values
        :type time: list
        :param val: values
        :type val: list
        :param err: errors, defaults to None
        :type err: list, optional
        :param ofac: over-factorization, defaults to 1
        :type ofac: int, optional
        :param normalize: normalization of periodogram, defaults to True
        :type normalize: bool, optional
        :param win_func: return window function, defaults to False
        :type win_func: bool, optional
        :return: frequency, power spectrum and phase of periodogram; if win_func is True, the frequency, power spectrum and phase of the window function are returned as well
        :rtype: tuple
        """

        # if not provided, set uncertainties to unity
        if err is None:
            err = np.ones(len(time))

        # frequencies, power spectrum and phases of data
        freq, ps, phi = _gls(time, val, err, ofac, normalize)

        # spectral window function
        if win_func:

            # central frequency
            freq_c = freq[int(len(freq)/2)]

            # sinusoid with unit amplitude at central frequency
            win_val = np.sin(2*np.pi*freq_c*time)
            win_err = err

            # power spectrum of window function
            win_freq, win_ps, _ = _gls(time, win_val, win_err, ofac)

            # recenter window function frequencies
            win_freq -= freq_c

            # area of window function
            win_dfreq = np.mean(win_freq[1:]-win_freq[:-1])
            win_area = np.sum(win_ps)*win_dfreq

        # return periodogram parameters
        if win_func:
            return freq, ps, phi, win_freq, win_ps, win_area
        else:
            return freq, ps, phi

def _gls(time, val, err, ofac, normalize=True):

    # time span and steps
    Time = time[-1] - time[0]
    dtime = time[1:] - time[:-1]

    # linear and angular frequencies
    dfreq = 1/(Time*ofac)
    freq = np.arange(1/Time, 1/(2*np.median(dtime)), dfreq)
    omega = 2*np.pi*freq

    # weights
    W = np.sum(1/err**2)
    w = 1/(W*err**2)

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
        #Y = _weight_sum1(w, val)
        #C = _weight_sum1(w, cosarg)
        #S = _weight_sum1(w, sinarg)
        Y = np.sum(w*val)
        C = np.sum(w*cosarg)
        S = np.sum(w*sinarg)

        # weighted sums of cross-terms
        #YYhat = _weight_sum2(w, val, val)
        #YChat = _weight_sum2(w, val, cosarg)
        #YShat = _weight_sum2(w, val, sinarg)
        #CChat = _weight_sum2(w, cosarg, cosarg)
        #SShat = _weight_sum2(w, sinarg, sinarg)
        #CShat = _weight_sum2(w, cosarg, sinarg)
        YYhat = np.sum(w*val*val)
        YChat = np.sum(w*val*cosarg)
        YShat = np.sum(w*val*sinarg)
        CChat = np.sum(w*cosarg*cosarg)
        SShat = np.sum(w*sinarg*sinarg)
        CShat = np.sum(w*cosarg*sinarg)

        # differences of sums
        YY = YYhat - Y * Y
        YC = YChat - Y * C
        YS = YShat - Y * S
        CC = CChat - C * C
        SS = SShat - S * S
        CS = CShat - C * S

        # normalization
        D = CC * SS - CS ** 2

        # amplitudes
        a = (YC * SS - YS * CS) / D
        b = (YS * CC - YC * CS) / D

        # power spectrum
        if normalize == True:
            ps[i] = (SS*YC**2 + CC*YS**2 - 2*CS*YC*YS) / (YY*D)
        if normalize == False:
            ps[i] = a**2 + b**2

        # phases
        phi[i] = np.arctan2(a, b)

    # return frequencies, power spectrum and phases
    return freq, ps, phi

#@njit()
def _weight_sum1(w, arr):
    sumw = 0
    for i in range(len(w)):
        sumw += w[i] * arr[i]
    return sumw

#@njit()
def _weight_sum2(w, arr1, arr2):
    sumw = 0
    for i in range(len(w)):
        sumw += w[i] * arr1[i] * arr2[i]
    return sumw