import numpy as np

class recovery_test:

    def recovery_test(self, P_inj:list, K_inj:list, p_inj:list=None, P_err:list=None, ofac:int=3, fap:float=0.01, N_max:int=10) -> list:
        """Recovery test (used by injection_recovery() function).

        :param P_inj: injected periods
        :type P_inj: list
        :param K_inj: injected RV semi-amplitudes
        :type K_inj: list
        :param p_inj: injected phases, defaults to None
        :type p_inj: list
        :param P_err: allowed period errors to count as recovery (set to 10% of injected periods if not provided), defaults to None
        :type P_err: list, optional
        :param ofac: over-factorization of periodogram, defaults to 3
        :type ofac: int, optional
        :param fap: false-alarm probability level, defaults to 0.01
        :type fap: float, optional
        :param N_max: maximum number of fitted Keplerians, defaults to 10
        :type N_max: int, optional
        :return: ratios between recovered and injected RV semi-amplitudes for recovered Keplerians, NaNs otherwise
        :rtype: list
        """

        # read data
        time_val, = [self.arve.data.time[key] for key in ["time_val"]]
        vrad_val, = [self.arve.data.vrad[key] for key in ["vrad_val"]]

        # nr. of injected Keplerians
        Nkep_inj = len(P_inj)

        # injected phase
        if p_inj is None:
            p_inj = np.random.uniform(-np.pi, np.pi, Nkep_inj)

        # period error
        if P_err is None:
            P_err = P_inj*0.1
        
        # copy RV values
        vrad_val_tmp = np.copy(vrad_val)

        # add injected Keplerians to RV values
        for i in range(Nkep_inj):
            vrad_val_tmp += K_inj[i]*np.sin(2*np.pi/P_inj[i]*time_val + p_inj[i])

        # set RV values with Keplerians
        self.arve.data.vrad["vrad_val"] = vrad_val_tmp

        # empty array for recovery results
        recovery_result = np.zeros(Nkep_inj)

        # try to fit Keplerians
        try:

            # fit Keplerians
            self.fit_keplerians(ofac=ofac, fap=fap, N_max=N_max)
            keplerians = self.keplerians
            Nkep_fit = len(keplerians)

            # get periods and RV semi-amplitudes of fitted Keplerians
            P_val = keplerians["P_val"].to_numpy()
            K_val = keplerians["K_val"].to_numpy()

            # loop injected Keplerians
            for i in range(Nkep_inj):

                # check if injected period is among fitted periods (within the allowed error)
                P_bound = np.vstack([P_val   -P_err  [i]    ,  P_val   +P_err  [i]]).T
                P_crit  = np.array([(P_inj[i]>P_bound[k,0]) & (P_inj[i]<P_bound[k,1]) for k in range(Nkep_fit)])

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

            recovery_result = np.zeros(Nkep_inj)*np.nan
        
        # re-set RV values without Keplerians
        self.arve.data.vrad["vrad_val"] = vrad_val

        # return recovery result
        return recovery_result