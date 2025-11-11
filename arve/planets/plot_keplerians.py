import matplotlib.pyplot as     plt
from   matplotlib.ticker import AutoMinorLocator
import numpy             as     np

class plot_keplerians:

    def plot_keplerians(
        self,
        figsize : tuple | None = None,
        ) -> plt.Figure:
        """Plot Keplerians.

        Parameters
        ----------
        figsize : tuple | None, optional
            figure size, by default None

        Returns
        -------
        plt.Figure
            figure with fitted Keplerians and their periodograms
        """

        # read data
        time_val,          = [self.arve.data.time[key] for key in ["time_val"]            ]
        vrad_val, vrad_err = [self.arve.data.vrad[key] for key in ["vrad_val", "vrad_err"]]

        # periodograms
        freq, power_gls_arr, power_fap_arr = [self.periodograms[key] for key in ["freq", "power_gls", "power_fap"]]

        # read Keplerians
        keplerians   = self.keplerians
        N_kep        = len(keplerians)
        cols_val     = ["val" in key for key in keplerians.keys()]
        cols_err     = ["err" in key for key in keplerians.keys()]
        para_val_arr = keplerians[keplerians.keys()[cols_val]].to_numpy()
        para_err_arr = keplerians[keplerians.keys()[cols_err]].to_numpy()

        # copy RV values
        vrad_val_tmp  = np.copy(vrad_val)

        # initiate figure
        if figsize is None:
            figsize = (15,(N_kep+1)*3)
        fig, axs = plt.subplots(N_kep+1, 2, figsize=figsize)

        # loop Keplerians
        for i in range(N_kep+1):            

            # plot phase-folded RV time series with fitted Keplerian
            if N_kep == 0:
                ax = axs[0]
            else:
                ax = axs[i,0]
            if i < N_kep:
                P = para_val_arr[i,0]
                ax.errorbar((time_val%P)/P, vrad_val_tmp, vrad_err, fmt=".", color="k", ecolor="r")
                vrad_val_mod_calc = self.arve.functions.keplerian(time_val, *para_val_arr[i])
                vrad_val_tmp     -= vrad_val_mod_calc
                time_val_plot     = np.linspace(0,P,int(1e4))
                vrad_val_mod_plot = self.arve.functions.keplerian(time_val_plot, *para_val_arr[i])
                ax.plot(time_val_plot/P, vrad_val_mod_plot, "-b", zorder=100)
            else:
                ax.errorbar(time_val, vrad_val_tmp, vrad_err, fmt=".", color="k", ecolor="r")
            
            # time series labels
            if i < N_kep-1:
                ax.set_xticklabels([])
                ax.set_ylabel("RV [km/s]")
            if i == N_kep-1:
                ax.set_xlabel("Phase")
                ax.set_ylabel("RV [km/s]")
            if i == N_kep:
                ax.set_xlabel("Time [d]")
                ax.set_ylabel("RV [km/s]")
            
            # time series axes
            ax.tick_params(axis="both", which="both", direction="in", top=True, right=True)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())

            # plot GLS periodogram
            if N_kep == 0:
                ax = axs[1]
            else:
                ax = axs[i,1]
            power_gls = power_gls_arr[i]
            power_fap = power_fap_arr[i]
            ax.semilogx(1/freq, power_gls, "-k")
            ax.axhline(power_fap, ls="--", c="k")
            if i < N_kep:
                ax.plot(1/freq[np.argmax(power_gls)], power_gls[np.argmax(power_gls)], "ob")
            
            # periodogram labels
            if i < N_kep-1:
                ax.set_xticklabels([])
                ax.set_ylabel("Power")
            if i == N_kep-1:
                ax.set_xlabel("Period [d]")
                ax.set_ylabel("Power")
            if i == N_kep:
                ax.set_xlabel("Period [d]")
                ax.set_ylabel("Power")
            
            # periodogram axes
            ax.set_ylim(0,ax.get_ylim()[1])
            ax.yaxis.set_label_position("right")
            ax.yaxis.tick_right()
            ax.tick_params(axis="both", which="both", direction="in", top=True, left=True)
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            
            # layout
            plt.tight_layout()

        # return figure
        return fig