import numpy as np

def add_vpsd_components(self, components:list=["Photon_noise", "Oscillations", "Granulation", "Supergranulation"]) -> None:
    """Add velocity power spectral density (VPSD) components.

    :param components: VPSD components, defaults to ["Photon_noise", "Oscillations", "Granulation", "Supergranulation"]
    :type components: list, optional
    :return: None
    :rtype: None
    """

    # read stellar parameters and VPSD
    Teff, logg, M, R = [self.stellar_parameters[key] for key in ["Teff", "logg", "M", "R"]]
    freq, vpsd, freq_avg, vpsd_avg = [self.vpsd[key] for key in ["freq", "vpsd", "freq_avg", "vpsd_avg"]]

    # photon noise
    if "Photon_noise" in components:

        # coefficients
        c0 = np.min(vpsd_avg) # amplitude

        # specifications
        name     = "Photon_noise" # name
        type     = "Constant"     # function type
        coef_val = [c0]           # coefficient values
        coef_err = [0]            # coefficient errors
        vary     = [True]         # vary when fitted

        # check noise is non-negative
        if c0 >= 0:

            # save VPSD component
            self.vpsd_components[name] = {"type": type, "coef_val": coef_val, "coef_err": coef_err, "vary": vary}

    # oscillations
    if "Oscillations" in components:

        # coefficients
        c2 = M/(R**2*np.sqrt(Teff/5777))*3.05e-3 * 60*60*24 # central frequency
        c1 = M**(1/2)*R**(-3/2)*134.9*1e-6*2 * 60*60*24     # peak width
        c0 = vpsd_avg[np.argmin(np.abs(freq_avg-c2))]       # amplitude

        # specifications
        name     = "Oscillations"     # name
        type     = "Lorentz"          # function type
        coef_val = [c0  , c1  , c2  ] # coefficient values
        coef_err = [0   , 0   , 0   ] # coefficient errors
        vary     = [True, True, True] # vary when fitted

        # check central frequency is within frequency range
        if (c2 > min(freq)) & (c2 < max(freq)):

            # save VPSD component
            self.vpsd_components[name] = {"type": type, "coef_val": coef_val, "coef_err": coef_err, "vary": vary}

    # granulation
    if "Granulation" in components:

        # coefficients
        c2 = 2                                                    # exponent
        c1 = 1/24*(10**logg/10**(4.4))**(7/9)*(Teff/5777)**(23/9) # characteristic timescale
        c0 = vpsd_avg[np.argmin(np.abs(freq_avg-1/c1))]           # amplitude

        # specifications
        name     = "Granulation"       # name
        type     = "Harvey"            # function type
        coef_val = [c0  , c1  , c2   ] # coefficient values
        coef_err = [0   , 0   , 0    ] # coefficient errors
        vary     = [True, True, False] # vary when fitted

        # check characteristic frequency is within frequency range
        if (1/c1 > min(freq)) & (1/c1 < max(freq)):

            # save VPSD component
            self.vpsd_components[name] = {"type": type, "coef_val": coef_val, "coef_err": coef_err, "vary": vary}

    # supergranulation
    if "Supergranulation" in components:

        # coefficients
        c2 = 2                                                     # exponent
        c1 = 10/24*(10**logg/10**(4.4))**(7/9)*(Teff/5777)**(23/9) # characteristic timescale
        c0 = vpsd_avg[np.argmin(np.abs(freq_avg-1/c1))]            # amplitude

        # specifications
        name     = "Supergranulation"  # name
        type     = "Harvey"            # function type
        coef_val = [c0  , c1  , c2   ] # coefficient values
        coef_err = [0   , 0   , 0    ] # coefficient errors
        vary     = [True, True, False] # vary when fitted

        # check characteristic frequency is within frequency range
        if (1/c1 > min(freq)) & (1/c1 < max(freq)):

            # save VPSD component
            self.vpsd_components[name] = {"type": type, "coef_val": coef_val, "coef_err": coef_err, "vary": vary}
    
    return None