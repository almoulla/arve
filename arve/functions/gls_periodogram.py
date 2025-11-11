import numpy as np

class gls_periodogram:

    def gls_periodogram(
        self,
        time_val  : np.ndarray               ,
        data_val  : np.ndarray               ,
        data_err  : np.ndarray | None = None ,
        oversamp  : float             = 1    ,
        normalize : bool              = True ,
        win_func  : bool              = False,
        ) -> tuple[np.ndarray | float]:
        """Generalized Lomb-Scargle (GLS) periodogram.

        Parameters
        ----------
        time_val : np.ndarray
            time values
        data_val : np.ndarray
            data values
        data_err : np.ndarray | None, optional
            data errors, by default None
        oversamp : float, optional
            oversamling factor of the periodogram frequency grid, by default 1
        normalize : bool, optional
            normalization of the periodogram, by default True
        win_func : bool, optional
            return window function, by default False

        Returns
        -------
        tuple[np.ndarray | float]
            frequency, power spectrum and phase of the periodogram; if win_func is True, the frequency, power spectrum and area of the window function are returned as well
        """

        # if not provided, set uncertainties to unity
        if data_err is None:
            data_err = np.ones(len(time_val))

        # frequencies, power spectrum and phases of data
        freq, ps, phi = _gls(time_val, data_val, data_err, oversamp, normalize)

        # spectral window function
        if win_func:

            # central frequency
            freq_c = freq[int(len(freq)/2)]

            # sinusoid with unit amplitude at central frequency
            win_val = np.sin(2*np.pi*freq_c*time_val)
            win_err = data_err

            # power spectrum of window function
            win_freq, win_ps, _ = _gls(time_val, win_val, win_err, oversamp)

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

def _gls(
    time_val  : np.ndarray ,
    data_val  : np.ndarray ,
    data_err  : np.ndarray ,
    oversamp  : float      ,
    normalize : bool = True
    ) -> tuple[np.ndarray]:
    """Generalized Lomb-Scargle (GLS) periodogram.

    Parameters
    ----------
    time_val : np.ndarray
        time values
    data_val : np.ndarray
        data values
    data_err : np.ndarray
        data errors
    oversamp : float
        oversamling factor of the periodogram frequency grid
    normalize : bool, optional
        normalization of the periodogram, by default True

    Returns
    -------
    tuple[np.ndarray]
        frequency, power spectrum and phase of the periodogram
    """

    # time parameters
    time_span = time_val[-1] - time_val[0]
    time_step = time_val[1:] - time_val[:-1]

    # frequency parameters
    freq_min  = 1/time_span
    freq_max  = 1/(2*np.median(time_step))
    freq_step = 1/(time_span*oversamp)
    freq      = np.arange(freq_min, freq_max, freq_step)
    omega     = 2*np.pi*freq

    # weights
    w  = 1/(data_err**2)
    w /= np.sum(w)

    # empty arrays for powers and phases
    N_freq = len(freq)
    ps     = np.empty(N_freq)
    phi    = np.empty(N_freq)

    # loop through frequencies
    for i in range(N_freq):

        # trigonometric terms
        arg    = omega[i]*time_val
        cosarg = np.cos(arg)
        sinarg = np.sin(arg)

        # weighted sums
        Y = np.sum(w*data_val)
        C = np.sum(w*cosarg)
        S = np.sum(w*sinarg)

        # weighted sums of cross-terms
        YYhat = np.sum(w*data_val*data_val)
        YChat = np.sum(w*data_val*cosarg)
        YShat = np.sum(w*data_val*sinarg)
        CChat = np.sum(w*cosarg*cosarg)
        SShat = np.sum(w*sinarg*sinarg)
        CShat = np.sum(w*cosarg*sinarg)

        # differences of sums
        YY = YYhat - Y*Y
        YC = YChat - Y*C
        YS = YShat - Y*S
        CC = CChat - C*C
        SS = SShat - S*S
        CS = CShat - C*S

        # normalization
        D = CC*SS - CS**2

        # amplitudes
        a = (YC*SS - YS*CS) / D
        b = (YS*CC - YC*CS) / D

        # power spectrum
        if normalize == True:
            ps[i] = (SS*YC**2 + CC*YS**2 - 2*CS*YC*YS) / (YY*D)
        if normalize == False:
            ps[i] = a**2 + b**2

        # phases
        phi[i] = np.arctan2(a, b)

    # return frequencies, power spectrum and phases
    return freq, ps, phi