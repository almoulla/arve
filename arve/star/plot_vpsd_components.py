import matplotlib.pyplot as plt
import numpy             as np

class plot_vpsd_components:

    def plot_vpsd_components(
        self,
        fig        : plt.Figure | None = None                                                        ,
        figsize    : tuple             = (10,10)                                                     ,
        colors     : list[str]         = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"],
        add_legend : bool              = True                                                        ,
        total_only : bool              = False                                                       ,
        ) -> plt.Figure:
        """Plot velocity power spectral density (VPSD) components.

        Parameters
        ----------
        fig : plt.Figure | None, optional
            figure on which to plot (if None, a new figure is created), by default None
        figsize : tuple, optional
            figure size, by default (10,10)
        colors : list[str], optional
            line colors for components, by default ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]
        add_legend : bool, optional
            add legend about line colors, by default True
        total_only : bool, optional
            plot only total component (sum of all components), by default False

        Returns
        -------
        plt.Figure
            figure with the VPSD, its log-average and its modeled components
        """

        # figure
        if fig is None:
            fig = plt.figure(figsize=figsize)

        # check VPSD is computed
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
        for i, comp_name in enumerate(self.vpsd_components.keys()):

            # component dictionary
            comp_dict = self.vpsd_components[comp_name]

            # function type and coefficients
            func_type = comp_dict["func_type"]
            coef_val  = comp_dict["coef_val"]

            # type Lorentz
            if func_type == "lorentz":
                
                # unpack coefficients
                c0 = coef_val[0]
                c1 = coef_val[1]
                c2 = coef_val[2]

                # compute component
                vpsd_comp = c0*c1**2/(c1**2+(freq-c2)**2)

            # type Harvey
            if func_type == "harvey":

                # unpack coefficients
                c0 = coef_val[0]
                c1 = coef_val[1]
                c2 = coef_val[2]

                # compute component
                vpsd_comp = c0/(1+(c1*freq)**c2)

            # type constant
            if func_type == "constant":

                # unpack coefficients
                c0 = coef_val[0]

                # compute component
                vpsd_comp = c0*np.ones_like(freq)

            # add component to sum
            vpsd_tot += vpsd_comp

            # plot component
            if total_only == False:
                plt.loglog(freq, vpsd_comp, ls="--", color=colors[i], label=comp_name.replace("_"," ").capitalize())

        # plot component sum
        if total_only == False:
            plt.loglog(freq, vpsd_tot, ls="-", c="k", label="Total")
        else:
            plt.loglog(freq, vpsd_tot, ls="-", c="k")
            plt.text(freq[0]*1.25, vpsd_tot[0]*0.85, self.stellar_parameters["sptype"], ha="left", va="top")

        # plot limits
        plt.xlim(min(freq), max(freq))

        # plot labels
        plt.xlabel("$f$ [d$^{-1}$]")
        plt.ylabel("VPSD [(km/s)$^{2}$ / d$^{-1}$]")

        # plot legend
        if (total_only == False) & (add_legend == True):
            leg = plt.legend(loc="lower left", framealpha=1)
            leg.set_zorder(101)

        # plot axes
        plt.gca().tick_params(axis="both", which="both", direction="in", top=True, right=True)

        # plot layout
        plt.tight_layout()

        # return figure
        return fig