import numpy as np

class detection_test:

    def detection_test(self, P_inj:float, K_inj:float, P_err:float=None, ofac:int=3, fap:float=0.01) -> float:
        """Detection test (used by injection_recovery() function).

        :param P_inj: injected period
        :type P_inj: float
        :param K_inj: injected RV semi-amplitude
        :type K_inj: float
        :param P_err: allowed period error to count as detection (set to 10% of injected period if not provided), defaults to None
        :type P_err: float, optional
        :param ofac: over-factorization of periodogram, defaults to 3
        :type ofac: int, optional
        :param fap: false-alarm probability level, defaults to 0.01
        :type fap: float, optional
        :return: ratio between recovered and injected RV semi-amplitude if detected, NaN otherwise
        :rtype: float
        """

        # read data
        time_val, = [self.arve.data.time[key] for key in ["time_val"]]
        vrad_val, = [self.arve.data.vrad[key] for key in ["vrad_val"]]

        # period error
        if P_err is None:
            P_err = P_inj*0.1

        try:

            # copy RV values
            vrad_val_tmp = np.copy(vrad_val)

            # add injected Keplerian to RV values
            vrad_val_tmp += K_inj*np.sin(2*np.pi/P_inj*time_val)

            # fit Keplerians
            self.arve.data.vrad["vrad_val"] = vrad_val_tmp
            self.fit_keplerians(ofac=ofac, fap=fap, P_err=P_err)
            keplerians = self.keplerians
            Nkepl = len(keplerians)

            # get periods and RV semi-amplitudes of fitted Keplerians
            P_val    = keplerians["P_val"].to_numpy()
            K_val    = keplerians["K_val"].to_numpy()

            # check if injected period is among fitted periods (within the allowed error)
            P_bound  = np.vstack([P_val-P_err, P_val+P_err]).T
            P_crit   = np.array([(P_inj>P_bound[k,0]) & (P_inj<P_bound[k,1]) for k in range(Nkepl)])

            # if period criterion is not satisfied, return NaN
            if np.sum(P_crit) == 0:
                detection_result = np.nan
            # if period criterion is satisfied, return ratio between recovered and injected RV semi-amplitude of the Keplerian with the closest period
            else:
                idx_all = np.where(P_crit)[0]
                idx_min = np.argmin(np.abs(K_val[idx_all]-K_inj))
                K_rec   = K_val[idx_all[idx_min]]
                detection_result = K_rec/K_inj
        
        except:

            detection_result = np.nan
        
        # re-set RV values without Keplerian
        self.arve.data.vrad["vrad_val"] = vrad_val

        # return detection result
        return detection_result