import numpy as np

class compute_vpsd:

    def compute_vpsd(
        self,
        N_bin : int | None = None,
        ) -> None:
        """Compute velocity power spectral density (VPSD).

        Parameters
        ----------
        N_bin : int | None, optional
            number of logarithmically equidistant frequency bins for the averaged VPSD (if None, each order of magnitude in frequency is divided into 10 bins), by default None

        Returns
        -------
        None
            None
        """

        # read RV time series
        time_val,          = [self.arve.data.time[key] for key in ["time_val"]]
        vrad_val, vrad_err = [self.arve.data.vrad[key] for key in ["vrad_val", "vrad_err"]]

        # compute velocity power spectrum
        freq, vps, phi, win_freq, win_vps, win_area = self.arve.functions.gls_periodogram(time_val=time_val, data_val=vrad_val, data_err=vrad_err, normalize=False, win_func=True)

        # compute log-average VPS
        if N_bin is None:
            N_bin = int((np.log10(freq[-1])-np.log10(freq[0]))*10)
        freq_bin = 10 ** (np.linspace(np.log10(freq[0]), np.log10(freq[-1]), N_bin+1))
        freq_avg = (freq_bin[1:] + freq_bin[:-1]) / 2
        vps_avg  = np.empty(freq_avg.size)
        for i in range(freq_avg.size):
            vps_bin = vps[(freq > freq_bin[i]) & (freq < freq_bin[i + 1])]
            if len(vps_bin) == 0:
                vps_avg[i] = np.nan
            else:
                vps_avg[i] = np.mean(vps_bin)

        # delete NaN values
        i_delete = np.isnan(vps_avg)
        freq_avg = np.delete(freq_avg, np.where(i_delete))
        vps_avg  = np.delete(vps_avg , np.where(i_delete))

        # compute VPSD
        vpsd     = vps    /win_area
        vpsd_avg = vps_avg/win_area

        # save VPSD
        self.vpsd = {"freq": freq, "vps": vps, "vpsd": vpsd, "phase": phi, "freq_avg": freq_avg, "vps_avg": vps_avg, "vpsd_avg": vpsd_avg}

        return None