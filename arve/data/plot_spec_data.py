import matplotlib.pyplot as     plt
from   matplotlib.ticker import AutoMinorLocator
import numpy             as     np

class plot_spec_data:

    def plot_spec_data(
        self,
        figsize          : tuple              = (20,10),
        xlim             : list[float] | None = None   ,
        ylim             : list[float] | None = None   ,
        orders           : list[int]   | None = None   ,
        include_obs_data : bool               = True   ,
        include_aux_data : bool               = True   ,
        plot_spec        : bool               = True   ,
        plot_tell        : bool               = True   ,
        plot_band        : bool               = True   ,
        plot_excl        : bool               = True   ,
        plot_mask        : bool               = True   ,
        add_legend       : bool               = True   ,
        annotate         : bool               = False  ,
        ) -> plt.Figure:
        """Plot spectral data.

        Parameters
        ----------
        figsize : tuple, optional
            figure size, by default (20,10)
        xlim : list[float] | None, optional
            x limits, by default None
        ylim : list[float] | None, optional
            y limits, by default None
        orders : list[int] | None, optional
            indices of orders to plot, by default None
        include_obs_data : bool, optional
            include observed data, by default True
        include_aux_data : bool, optional
            include auxiliary data, by default True
        plot_spec : bool, optional
            plot stellar spectrum, by default True
        plot_tell : bool, optional
            plot telluric spectrum, by default True
        plot_band : bool, optional
            plot telluric bands, by default True
        plot_excl : bool, optional
            plot excluded wavelength regions, by default True
        plot_mask : bool, optional
            plot stellar line mask, by default True
        add_legend : bool, optional
            add legend about line colors, by default True
        annotate : bool, optional
            annonate spectral lines with atomic information, by default False

        Returns
        -------
        plt.Figure
            figure with the spectral data
        """
        

        # read data
        wave_val, flux_val, flux_err = [self.spec    [key] for key in ["wave_val", "flux_val", "flux_err"    ]]
        spec, tell, band, excl, mask = [self.aux_data[key] for key in ["spec", "tell", "band", "excl", "mask"]]
        if flux_val is not None:
            flux_val = flux_val[0]
        else:
            flux_val = np.zeros_like(wave_val)
        if flux_err is not None:
            flux_err = flux_err[0]
        else:
            flux_err = np.zeros_like(wave_val)
        N_ord = wave_val.shape[0]
        
        # orders to plot
        if orders is None:
            orders = np.arange(0, N_ord)
        
        # xlimits
        if xlim is None:
            xlim = [np.nanmin(wave_val[orders]), np.nanmax(wave_val[orders])]

        # redefining orders for efficient plotting
        if len(orders) == N_ord:
            ord_min = orders[0]
            while wave_val[ord_min][-1]<=xlim[ 0]:
                ord_min += 1
            ord_max = ord_min
            while wave_val[ord_max][ 0]<=xlim[-1]:
                ord_max += 1
                if ord_max == N_ord:
                    break
            orders = np.arange(ord_min, ord_max)

        # initate figure
        fig = plt.figure(figsize=figsize)
        if include_obs_data & include_aux_data:
            ax1  = plt.gca()
            ax2  = plt.gca().twinx()
            axs  = [ax1, ax2]
            unit = ["[observed units]", "[normalized]"]
        elif include_obs_data:
            ax1  = plt.gca()
            axs  = [ax1]
            unit = ["[observed units]"]
        elif include_aux_data:
            ax2  = plt.gca()
            axs  = [ax2]
            unit = ["[normalized]"]
        
        # include observed data
        if include_obs_data:

            # loop orders
            for o in orders:

                # index of points within x limits
                idx_xlim = (wave_val[o]>=xlim[0]) & (wave_val[o]<=xlim[1])
            
                # plot observed stellar spectrum
                ax1.errorbar(wave_val[o][idx_xlim], flux_val[o][idx_xlim], flux_err[o][idx_xlim], capsize=1, zorder=3, color="k", ecolor="r")

        # include auxiliary data
        if include_aux_data:

            # loop orders
            for o in orders:

                # index of points within x limits
                idx_xlim = (spec["wave"][o]>=xlim[0]) & (spec["wave"][o]<=xlim[1])

                # plot stellar spectrum
                if plot_spec:
                    ax2.plot(spec["wave"][o][idx_xlim], spec["flux"][o][idx_xlim], color="k", zorder=3)

                # plot telluric spectrum
                if plot_tell:
                    ax2.plot(tell["wave"][o][idx_xlim], tell["flux"][o][idx_xlim], color="g", zorder=2)

                # plot telluric bands
                if plot_band:
                    for i in range(len(band[o])):
                        if (band[o]["wave_u"][i]>=xlim[0]) & (band[o]["wave_l"][i]<=xlim[1]):
                            ax2.axvspan(band[o]["wave_l"][i], band[o]["wave_u"][i], alpha=0.25, color="g", zorder=0)
                
                # plot excluded wavelength regions
                if plot_excl:
                    for i in range(len(excl[o])):
                        if (excl[o]["wave_u"][i]>=xlim[0]) & (excl[o]["wave_l"][i]<=xlim[1]):
                            ax2.axvspan(excl[o]["wave_l"][i], excl[o]["wave_u"][i], alpha=0.25, color="r", zorder=0)
                
                # plot stellar line mask
                if plot_mask:
                    for i in range(len(mask[o])):
                        if (mask[o]["wave_u"][i]>=xlim[0]) & (mask[o]["wave_l"][i]<=xlim[1]):
                            if    mask[o]["crit_tell"][i] &  mask[o]["crit_excl"][i]:
                                color = "k"
                            elif ~mask[o]["crit_tell"][i] &  mask[o]["crit_excl"][i]:
                                color = "g"
                            elif  mask[o]["crit_tell"][i] & ~mask[o]["crit_excl"][i]:
                                color = "r"
                            elif ~mask[o]["crit_tell"][i] & ~mask[o]["crit_excl"][i]:
                                color = "y"
                            ax2.vlines (mask[o]["wave"  ][i], mask[o]["flux"  ][i], 1, ls="--", color=color, zorder=1)
                            ax2.axvspan(mask[o]["wave_l"][i], mask[o]["wave_u"][i], alpha=0.25, color="k"  , zorder=0)
                            if annotate:
                                annotation = mask[o]["species"][i]
                                annotation = annotation.replace(" ", "")
                                for number, roman in {"1": "I", "2": "II", "3": "III", "4": "IV", "5": "V", "6": "VI", "7": "VII", "8": "VIII", "9": "IX"}.items():
                                    annotation = annotation.replace(number, roman)
                                ax2.annotate(annotation, (mask[o]["wave"][i], 1.15), ha="center", va="center", fontfamily="monospace", fontsize=7)
            
            # plot legend
            if add_legend:
                ax2.plot([], [], ls="--", color="k", label="Unaffected")
                ax2.plot([], [], ls="--", color="g", label="Overlapping with telluric region")
                ax2.plot([], [], ls="--", color="r", label="Overlapping with excluded region")
                ax2.plot([], [], ls="--", color="y", label="Overlapping with both telluric and excluded region")
                ax2.legend(title="Stellar lines", loc="lower left", framealpha=1).set_zorder(100)

        # plot labels
        for i in range(len(axs)):
            axs[i].set_xlabel("$\lambda$ [Å]")
            axs[i].set_ylabel("Flux " + unit[i])

        # plot axes
        for ax in axs:
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            if len(axs) == 1:
                ax.tick_params(axis="both", which="both", direction="in", top=True, right=True)
            else:
                ax.tick_params(axis="both", which="both", direction="in", top=True)

        # plot layout
        plt.tight_layout()

        # return figure
        return fig