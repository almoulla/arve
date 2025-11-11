from   astropy.timeseries import LombScargle
import numpy              as     np
import pandas             as     pd
from   scipy.optimize     import curve_fit

class fit_keplerians:

    def fit_keplerians(
        self,
        oversamp : float = 3   ,
        fap      : float = 0.01,
        P_lim    : float = 0.1 ,
        N_max    : int   = 10  ,
        ) -> None:
        """Fit Keplerians.

        Parameters
        ----------
        oversamp : float, optional
            oversamling factor of the periodogram frequency grid, by default 3
        fap : float, optional
            false-alarm probability (FAP) level, by default 0.01
        P_lim : float, optional
            allowed fractional period error for fitting bound, by default 0.1
        N_max : int, optional
            maximum number of fitted Keplerians, by default 10

        Returns
        -------
        None
            None
        """

        # read data
        time_val,          = [self.arve.data.time[key] for key in ["time_val"]            ]
        vrad_val, vrad_err = [self.arve.data.vrad[key] for key in ["vrad_val", "vrad_err"]]

        # copy RV values
        vrad_val_tmp = np.copy(vrad_val)

        # time parameters
        time_span = time_val[-1] - time_val[0]
        time_step = time_val[1:] - time_val[:-1]

        # frequency parameters
        freq_min  = 1/time_span
        freq_max  = 1/(2*np.median(time_step))
        freq_step = 1/(time_span*oversamp)
        freq      = np.arange(freq_min, freq_max, freq_step)

        # empty lists for parameter values and errors of fitted Keplerians
        para_val_arr = []
        para_err_arr = []

        # empty lists for periodogram power and FAP power
        power_gls_arr = []
        power_fap_arr = []

        # initial value for maximum periodogram power and FAP power
        power_max = 0
        power_fap = 0

        # fit Keplerians while the maximum periodogram power is larger than the FAP power or while the number of fitted Keplerians is below the maximum
        while (power_max >= power_fap) & (len(para_val_arr) < N_max):

            # Lomb-Scargle periodogram
            gls       = LombScargle(time_val, vrad_val_tmp, vrad_err)
            power_gls = gls.power(freq)
            power_max = np.nanmax(power_gls)
            power_fap = gls.false_alarm_level(fap)
            power_gls_arr.append(power_gls)
            power_fap_arr.append(power_fap)

            # if the maximum periodogram is larger than the FAP power, fit a Keplerian
            if power_max >= power_fap:

                # parameter guesses
                P_guess = 1/freq[np.argmax(power_gls)]
                K_guess = np.nanmax(vrad_val_tmp) - np.nanmedian(vrad_val_tmp)
                p_guess = 0
                C_guess = np.nanmedian(vrad_val_tmp)
                p0      = [P_guess, K_guess, p_guess, C_guess]

                # set period error if not provided
                P_err = P_guess*P_lim

                # parameter bounds
                P_bound = [P_guess-P_err, P_guess+P_err]
                K_bound = [0, np.inf]
                p_bound = [-np.pi, np.pi]
                C_bound = [-np.inf, np.inf]
                bounds  = [P_bound, K_bound, p_bound, C_bound]
                bounds  = np.transpose(bounds).tolist()

                # fit Keplerian, store parameter values and errors, and subtract Keplerian from RV values
                para_val, para_err = curve_fit(self.arve.functions.keplerian, time_val, vrad_val_tmp, sigma=vrad_err, p0=p0, bounds=bounds)
                para_err = np.diag(para_err)**(1/2)
                vrad_val_mod = self.arve.functions.keplerian(time_val, *para_val)
                para_val_arr.append(para_val)
                para_err_arr.append(para_err)
                vrad_val_tmp -= vrad_val_mod
        
        # Lomb-Scargle periodogram of residuals
        gls       = LombScargle(time_val, vrad_val_tmp, vrad_err)
        power_gls = gls.power(freq)
        power_max = np.nanmax(power_gls)
        power_fap = gls.false_alarm_level(fap)
        power_gls_arr.append(power_gls)
        power_fap_arr.append(power_fap)

        # convert periodograms to NumPy arrays
        power_gls_arr = np.array(power_gls_arr)
        power_fap_arr = np.array(power_fap_arr)

        # convert Keplerian parameters to NumPy arrays
        para_val_arr = np.array(para_val_arr)
        para_err_arr = np.array(para_err_arr)

        # store Keplerians as DataFrame
        keplerians = pd.DataFrame()
        if len(para_val_arr) > 0:
            keplerians["P_val"] = para_val_arr[:,0]
            keplerians["P_err"] = para_err_arr[:,0]
            keplerians["K_val"] = para_val_arr[:,1]
            keplerians["K_err"] = para_err_arr[:,1]
            keplerians["p_val"] = para_val_arr[:,2]
            keplerians["p_err"] = para_err_arr[:,2]
            keplerians["C_val"] = para_val_arr[:,3]
            keplerians["C_err"] = para_err_arr[:,3]

        # save periodograms and Keplerians
        self.periodograms = {
            "freq"     : freq,
            "power_gls": power_gls_arr,
            "power_fap": power_fap_arr,
            "fap"      : fap
        }
        self.keplerians = keplerians

        return None