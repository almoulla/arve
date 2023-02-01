"""
add vpsd components
"""

import numpy as np


def add_vpsd_components(star, components=["Photon_noise", "Oscillations", "Granulation", "Supergranulation"]):

    # read stellar parameters and VPSD
    Teff, logg, M, R = [star.stellar_parameters[var] for var in ["Teff", "logg", "M", "R"]]
    freq, vpsd, freq_avg, vpsd_avg = [star.vpsd[var] for var in ["freq", "vpsd", "freq_avg", "vpsd_avg"]]

    if "Photon_noise" in components:

        c0 = np.min(vpsd_avg) # amplitude

        name     = "Photon_noise"
        type     = "Constant"
        coef_val = [c0]
        coef_err = [0]
        vary     = [True]

        # save VPSD component
        star.vpsd_components[name] = {"type": type, "coef_val": coef_val, "coef_err": coef_err, "vary": vary}

    if "Oscillations" in components:

        c2 = M/(R**2*np.sqrt(Teff/5777))*3.05e-3 * 60*60*24 # center
        c1 = M**(1/2)*R**(-3/2)*134.9*1e-6*2 * 60*60*24     # width
        c0 = vpsd_avg[np.argmin(np.abs(freq_avg-c2))]       # amplitude

        name     = "Oscillations"
        type     = "Lorentz"
        coef_val = [c0, c1, c2]
        coef_err = [0,0,0]
        vary     = [True, True, True]

        # check within frequency range
        if (c2 > min(freq)) & (c2 < max(freq)):

            # save VPSD component
            star.vpsd_components[name] = {"type": type, "coef_val": coef_val, "coef_err": coef_err, "vary": vary}

    if "Granulation" in components:

        c2 = 2                                                        # exponent
        c1 = 1/24 * (10**logg/10**(4.4))**(7/9) * (Teff/5777)**(23/9) # timescale
        c0 = vpsd_avg[np.argmin(np.abs(freq_avg-1/c1))]               # amplitude

        name     = "Granulation"
        type     = "Harvey"
        coef_val = [c0, c1, c2]
        coef_err = [0,0,0]
        vary     = [True, True, False]

        # check within frequency range
        if (1/c1 > min(freq)) & (1/c1 < max(freq)):

            # save VPSD component
            star.vpsd_components[name] = {"type": type, "coef_val": coef_val, "coef_err": coef_err, "vary": vary}

    if "Supergranulation" in components:

        c2 = 2                                          # exponent
        c1 = c1*10                                      # timescale
        c0 = vpsd_avg[np.argmin(np.abs(freq_avg-1/c1))] # amplitude

        name     = "Supergranulation"
        type     = "Harvey"
        coef_val = [c0, c1, c2]
        coef_err = [0,0,0]
        vary     = [True, True, False]

        # check within frequency range
        if (1/c1 > min(freq)) & (1/c1 < max(freq)):

            # save VPSD component
            star.vpsd_components[name] = {"type": type, "coef_val": coef_val, "coef_err": coef_err, "vary": vary}