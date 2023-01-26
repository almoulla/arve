"""
gls periodogram
"""

import matplotlib.pyplot as plt
import numpy as np
from numba import njit


def gls_periodogram(functions, time, y, sig=0, ofac=1, normalize=True, win_func=False):

    # if not provided, set uncertainties to unity
    if type(sig) == int:
        sig = np.ones(len(time))

    # power spectrum and phases of data
    f, ps, phi = gls(time, y, sig, ofac, normalize)

    # spectral window function
    if win_func:

        # central frequency
        fc = f[int(len(f) / 2)]

        # sinusoid with unit amplitude at central frequency
        win_y = np.sin(2 * np.pi * fc * time)

        # power spectrum of window function
        win_f, win_ps, _ = gls(time, win_y, sig, ofac)

        # recenter window function frequencies
        win_f -= fc

        # area of window function
        win_df = np.mean(win_f[1:] - win_f[:-1])
        win_area = np.sum(win_ps) * win_df

    # return periodogram parameters
    if win_func:
        return f, ps, phi, win_f, win_ps, win_area
    else:
        return f, ps, phi


def gls(t, y, sig, ofac, normalize=True):

    # time span and steps
    T = t[-1] - t[0]
    dt = t[1:] - t[:-1]

    # linear and angular frequencies
    df = 1 / (T * ofac)
    f = np.arange(1 / T, 1 / (2 * np.median(dt)), df)
    omega = 2 * np.pi * f

    # weights
    W = sum(1 / sig ** 2)
    w = 1 / (W * sig ** 2)

    # empty arrays for powers and phases
    N = len(f)
    ps = np.empty(N)
    phi = np.empty(N)

    # loop through frequencies
    for i in range(N):

        # trigonometric terms
        arg = omega[i] * t
        cosarg = np.cos(arg)
        sinarg = np.sin(arg)

        # weighted sums
        Y = weight_sum1(w, y)
        C = weight_sum1(w, cosarg)
        S = weight_sum1(w, sinarg)

        # weighted sums of cross-terms
        YYhat = weight_sum2(w, y, y)
        YChat = weight_sum2(w, y, cosarg)
        YShat = weight_sum2(w, y, sinarg)
        CChat = weight_sum2(w, cosarg, cosarg)
        SShat = weight_sum2(w, sinarg, sinarg)
        CShat = weight_sum2(w, cosarg, sinarg)

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
            ps[i] = (SS * YC ** 2 + CC * YS ** 2 - 2 * CS * YC * YS) / (YY * D)
        if normalize == False:
            ps[i] = a ** 2 + b ** 2

        # phases
        phi[i] = np.arctan2(a, b)

    # return power spectrum and phases
    return f, ps, phi


@njit()
def weight_sum1(w, arr):
    sumw = 0
    for i in range(len(w)):
        sumw += w[i] * arr[i]
    return sumw


@njit()
def weight_sum2(w, arr1, arr2):
    sumw = 0
    for i in range(len(w)):
        sumw += w[i] * arr1[i] * arr2[i]
    return sumw