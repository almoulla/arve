import matplotlib.pyplot as plt
import numpy             as np

class plot_vpsd_components:

    def plot_vpsd_components(self, fig:plt.Figure=None, total_only:bool=False) -> plt.Figure:
        """Plot velocity power spectral density (VPSD) components.

        :param fig: figure on which to plot (if None, a new figure is created), defaults to None
        :type fig: plt.Figure, optional
        :param total_only: plot only total component (sum of all components), defaults to False
        :type total_only: bool, optional
        :return: figure with the VPSD, its log-average and its modeled components
        :rtype: plt.Figure
        """

        # figure
        if fig is None:
            fig = plt.figure()

        # check vpsd is computed
        if self.vpsd is not None:

            # read VPSD and units
            freq, vpsd, freq_avg, vpsd_avg = [self.vpsd[key] for key in ["freq", "vpsd", "freq_avg", "vpsd_avg"]]

            # plot VPSD and average VPSD
            plt.loglog(freq, vpsd, ls="-", c="k", alpha=0.5)
            plt.loglog(freq_avg, vpsd_avg, ls="None", marker="o", mec="k", mfc="None")
        
        else:

            # example frequency grid
            T  = 10
            dt = 1/(24*60*2)
            freq = np.arange(1/T, 1/(2*dt), 1/T)

        # empty array for sum of components
        vpsd_tot = np.zeros(len(freq))

        # loop components
        for comp in self.vpsd_components.keys():

            # component dictionary
            comp_dict = self.vpsd_components[comp]

            # type and coefficients
            type     = comp_dict["type"]
            coef_val = comp_dict["coef_val"]

            # type Lorentz
            if type == "Lorentz":
                
                # unpack coefficients
                c0 = coef_val[0]
                c1 = coef_val[1]
                c2 = coef_val[2]

                # compute component
                vpsd_comp = c0*c1**2/(c1**2+(freq-c2)**2)

                # plot component
                if total_only == False:
                    plt.loglog(freq, vpsd_comp, ls="--", label=comp)

            # type Harvey
            if type == "Harvey":

                # unpack coefficients
                c0 = coef_val[0]
                c1 = coef_val[1]
                c2 = coef_val[2]

                # compute component
                vpsd_comp = c0/(1+(c1*freq)**c2)

                # plot component
                if total_only == False:
                    plt.loglog(freq, vpsd_comp, ls="--", label=comp)

            # type Constant
            if type == "Constant":

                # unpack coefficients
                c0 = coef_val[0]

                # compute component
                vpsd_comp = c0

                # plot component
                if total_only == False:
                    plt.axhline(vpsd_comp, ls="--", label=comp)

            # add component to sum
            vpsd_tot += vpsd_comp

        # plot component sum
        if total_only == False:
            plt.loglog(freq, vpsd_tot, ls="-", c="k", label="Total")
        else:
            plt.loglog(freq, vpsd_tot, ls="-", c="k")
            plt.text(freq[0]*1.25, vpsd_tot[0]*1.1, self.stellar_parameters["sptype"], ha="left", va="bottom")

        # plot limits
        plt.xlim(min(freq), max(freq))

        # plot labels
        plt.xlabel("$f$ [d$^{-1}$]")
        plt.ylabel("VPSD [(km/s)$^{2}$ / d$^{-1}$]")

        # plot legend
        if total_only == False:
            leg = plt.legend(loc="lower left")
            leg.set_zorder(101)

        # plot axes
        plt.gca().tick_params(axis="both", which="both", direction="in", top=True, right=True)

        # plot layout
        plt.tight_layout()

        # return figure
        return fig