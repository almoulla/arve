import numpy as np

class recovery_test:

    def recovery_test(
        self,
        P_inj    : np.ndarray              ,
        K_inj    : np.ndarray              ,
        p_inj    : np.ndarray | None = None,
        P_err    : np.ndarray | None = None,
        oversamp : float             = 3   ,
        fap      : float             = 0.01,
        N_max    : int               = 10  ,
        ) -> np.ndarray:
        """Recovery test (used by injection_recovery() function).

        Parameters
        ----------
        P_inj : np.ndarray
            injected periods
        K_inj : np.ndarray
            injected RV semi-amplitudes
        p_inj : np.ndarray | None, optional
            injected phases, by default None
        P_err : np.ndarray | None, optional
            allowed differences between injected and fitted periods to count as recoveries (set to 10% of injected periods if not provided), by default None
        oversamp : float, optional
            oversamling factor of the periodogram frequency grid, by default 3
        fap : float, optional
            false-alarm probability (FAP) level, by default 0.01
        N_max : int, optional
            maximum number of fitted Keplerians, by default 10

        Returns
        -------
        np.ndarray
            ratios between recovered and injected RV semi-amplitudes for recovered Keplerians
        """

        # read data
        time_val, = [self.arve.data.time[key] for key in ["time_val"]]
        vrad_val, = [self.arve.data.vrad[key] for key in ["vrad_val"]]

        # nr. of injected Keplerians
        N_inj = len(P_inj)

        # injected phase
        if p_inj is None:
            p_inj = np.random.uniform(-np.pi, np.pi, N_inj)

        # period error
        if P_err is None:
            P_err = P_inj*0.1
        
        # copy RV values
        vrad_val_tmp = np.copy(vrad_val)

        # add injected Keplerians to RV values
        for i in range(N_inj):
            vrad_val_tmp += K_inj[i]*np.sin(2*np.pi/P_inj[i]*time_val + p_inj[i])

        # set RV values with Keplerians
        self.arve.data.vrad["vrad_val"] = vrad_val_tmp

        # empty array for recovery results
        recovery_result = np.zeros(N_inj)

        # try to fit Keplerians
        try:

            # fit Keplerians
            self.fit_keplerians(oversamp=oversamp, fap=fap, N_max=N_max)
            keplerians = self.keplerians
            N_fit = len(keplerians)

            # get periods and RV semi-amplitudes of fitted Keplerians
            P_val = keplerians["P_val"].to_numpy()
            K_val = keplerians["K_val"].to_numpy()

            # loop injected Keplerians
            for i in range(N_inj):

                # check if injected period is among fitted periods (within the allowed error)
                P_bound = np.vstack([P_val   -P_err  [i]    ,  P_val   +P_err  [i]]).T
                P_crit  = np.array([(P_inj[i]>P_bound[k,0]) & (P_inj[i]<P_bound[k,1]) for k in range(N_fit)])

                # if period criterion is not satisfied, return NaN
                if np.sum(P_crit) == 0:
                    recovery_result[i] = np.nan
                # if period criterion is satisfied, return ratio between recovered and injected RV semi-amplitude of the Keplerian with the closest period
                else:
                    idx_all = np.where(P_crit)[0]
                    idx_min = np.argmin(np.abs(K_val[idx_all]-K_inj[i]))
                    K_rec   = K_val[idx_all[idx_min]]
                    recovery_result[i] = K_rec/K_inj[i]

        # if unable to fit Keplerians, return NaNs
        except:

            recovery_result = np.zeros(N_inj)*np.nan
        
        # re-set RV values without Keplerians
        self.arve.data.vrad["vrad_val"] = vrad_val

        # return recovery result
        return recovery_result