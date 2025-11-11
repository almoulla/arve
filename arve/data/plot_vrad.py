import matplotlib.pyplot as     plt
from   matplotlib.ticker import AutoMinorLocator
import numpy             as     np

class plot_vrad:

    def plot_vrad(
        self,
        figsize : tuple = (20,10),
        ) -> plt.Figure:
        """Plot RV time series.

        Parameters
        ----------
        figsize : tuple, optional
            figure size, by default (20,10)

        Returns
        -------
        plt.Figure
            figure with the RV time series
        """

        # read data
        time_val,          = [self.time[key] for key in ["time_val"]]
        vrad_val, vrad_err = [self.vrad[key] for key in ["vrad_val", "vrad_err"]]

        # initate figure
        fig = plt.figure(figsize=figsize)

        # plot RV time series
        plt.errorbar(time_val, vrad_val-np.nanmedian(vrad_val), vrad_err, fmt=".", color="k", ecolor="r")

        # plot labels
        plt.xlabel("Time [d]")
        plt.ylabel("RV [km/s]")

        # plot axes
        plt.gca().tick_params(axis="both", which="both", direction="in", top=True, right=True)
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator())
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator())

        # plot layout
        plt.tight_layout()

        # return figure
        return fig