"""
plot vpsd components
"""

import matplotlib.pyplot as plt
import numpy as np


def plot_vpsd_components(star, label):

    # figure
    fig = plt.figure()

    # read VPSD and units
    freq, vpsd = [star.vpsd[label][var] for var in ["freq", "vpsd"]]
    time_unit, rv_unit = [star.arve.data.rv[label][var] for var in ["time_unit", "rv_unit"]]

    # log-average VPSD
    freq_bin = 10 ** (np.linspace(np.log10(freq[0]), np.log10(freq[-1]), 51))
    freq_avg = (freq_bin[1:] + freq_bin[:-1]) / 2
    vpsd_avg = np.empty(freq_avg.size)
    for i in range(freq_avg.size):
        vpsd_bin = vpsd[(freq > freq_bin[i]) & (freq < freq_bin[i + 1])]
        if len(vpsd_bin) == 0:
            vpsd_avg[i] = np.nan
        else:
            vpsd_avg[i] = np.mean(vpsd_bin)

    # delete NaN values
    i_delete = np.isnan(vpsd_avg)
    freq_avg = np.delete(freq_avg, np.where(i_delete))
    vpsd_avg = np.delete(vpsd_avg, np.where(i_delete))

    # plot VPSD and average VPSD
    plt.loglog(freq, vpsd, ls="-", c="k", alpha=0.5)
    plt.loglog(freq_avg, vpsd_avg, ls="None", marker="o", mec="k", mfc="None")

    # empty array for sum of components
    vpsd_tot = np.zeros(len(freq))

    # loop components
    for comp in star.vpsd_components[label].keys():

        # component dictionary
        comp_dict = star.vpsd_components[label][comp]

        # type and coefficients
        type = comp_dict["type"]
        coef = comp_dict["coef"]

        # component Constant
        if type == "Constant":

            # unpack coefficients
            c0 = coef[0]

            # compute component
            vpsd_comp = c0

            # plot component
            plt.axhline(vpsd_comp, ls="--", label=comp)

        # component Lorentz
        if type == "Lorentz":
            
            # unpack coefficients
            c0 = coef[0]
            c1 = coef[1]
            c2 = coef[2]

            # compute component
            vpsd_comp = c0*c1**2/(c1**2+(freq-c2)**2)

            # plot component
            plt.loglog(freq, vpsd_comp, ls="--", label=comp)

        # component Harvey
        if type == "Harvey":

            # unpack coefficients
            c0 = coef[0]
            c1 = coef[1]
            c2 = coef[2]

            # compute component
            vpsd_comp = c0/(1+(c1*freq)**c2)

            # plot component
            plt.loglog(freq, vpsd_comp, ls="--", label=comp)

        # add component to sum
        vpsd_tot += vpsd_comp

    # plot component sum
    plt.loglog(freq, vpsd_tot, ls="-", c="k", label="Total")

    # plot limits
    plt.xlim(freq[0], freq[-1])

    # plot labels
    plt.xlabel("$f$" + " " + f"[{time_unit}" + "$^{-1}$]")
    plt.ylabel("VPSD" + " " + f"[({rv_unit})" + "$^{2}$" + " / " + f"{time_unit}" + "$^{-1}$]")

    # plot legend
    leg = plt.legend(loc="lower left")
    leg.set_zorder(101)

    # plot layout
    plt.tight_layout()

    # return figure
    return fig