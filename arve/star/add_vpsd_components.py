import numpy as np

class add_vpsd_components:

    def add_vpsd_components(self, components:list=["Oscillations", "Granulation", "Supergranulation", "Noise"]) -> None:
        """Add velocity power spectral density (VPSD) components.

        :param components: VPSD components, defaults to ["Oscillations", "Granulation", "Supergranulation", "Noise"]
        :type components: list, optional
        :return: None
        :rtype: None
        """

        # read stellar parameters and VPSD
        Teff, logg, M, R = [self.stellar_parameters[key] for key in ["Teff", "logg", "M", "R"]]
        g = 10**logg
        if self.vpsd is not None:
            freq, vpsd, freq_avg, vpsd_avg = [self.vpsd[key] for key in ["freq", "vpsd", "freq_avg", "vpsd_avg"]]

        # stellar parameters for the Sun
        Teff_sun = 5770
        logg_sun = 4.4
        g_sun    = 10**logg_sun

        # constants
        microHz_to_Hz = 1e-6                                        # micro-Hertz        to Hertz
        milliHz_to_Hz = 1e-3                                        # milli-Hertz        to Hertz
        mps_to_kmps   = 1e-3                                        # m/s                to km/s
        sec_to_day    = 1/(60*60*24)                                # seconds            to days
        std_to_fwhm   = 2*np.sqrt(2*np.log(2))                      # standard deviation to full-width at half-maxmimum
        photo_to_velo = (Teff/Teff_sun)**(32/9) * (g/g_sun)**(-2/9) # photometry         to velocimetry scaling ratio

        # oscillations maximum (central) frequency
        nu_max_sun = 1e3/3.66
        nu_max     = nu_max_sun * (Teff/Teff_sun)**(-0.5) * (g/g_sun)

        # dictionary with VPSD components
        if self.vpsd_components is None:
            self.vpsd_components = {}

        # oscillations
        if "Oscillations" in components:

            # coefficients for the Sun
            c2_sun = nu_max_sun                           # central frequency
            c1_sun = 2.89e-1 * milliHz_to_Hz / sec_to_day # peak width
            c0_sun = 1.25e3 * mps_to_kmps**2 * sec_to_day # amplitude

            # coefficients
            c2 = nu_max                                                                         # central frequency
            c1 = c1_sun * (nu_max/nu_max_sun)**0.88                                             # peak width
            c0 = c0_sun * (nu_max/nu_max_sun)**(-0.609*2/0.855-0.88) * photo_to_velo**(2/0.855) # amplitude

            # coefficients
            #c2 = M/(R**2*np.sqrt(Teff/5777))*3.05e-3 * 60*60*24 # central frequency
            #c1 = M**(1/2)*R**(-3/2)*134.9*1e-6*2 * 60*60*24     # peak width
            #c0 = vpsd_avg[np.argmin(np.abs(freq_avg-c2))]       # amplitude

            # specifications
            name     = "Oscillations"     # name
            type     = "Lorentz"          # function type
            coef_val = [c0  , c1  , c2  ] # coefficient values
            coef_err = [0   , 0   , 0   ] # coefficient errors
            vary     = [True, True, True] # vary when fitted

            # check vpsd is computed
            if self.vpsd is not None:

                # check central frequency is within frequency range
                if (c2 > min(freq)) & (c2 < max(freq)):

                    # save VPSD component
                    self.vpsd_components[name] = {"type": type, "coef_val": coef_val, "coef_err": coef_err, "vary": vary}
            
            else:

                # save VPSD component
                self.vpsd_components[name] = {"type": type, "coef_val": coef_val, "coef_err": coef_err, "vary": vary}

        # granulation
        if "Granulation" in components:

            # coefficients for the Sun
            c2_sun = 2                                    # exponent
            c1_sun = 9.11e-1 / 24                         # characteristic timescale
            c0_sun = 4.11e2 * mps_to_kmps**2 * sec_to_day # amplitude

            # coefficients
            c2 = c2_sun                                                           # exponent
            c1 = c1_sun * (nu_max/nu_max_sun)**(-0.97)                            # characteristic timescale
            c0 = c0_sun * (nu_max/nu_max_sun)**(-0.609*2-0.97) * photo_to_velo**2 # amplitude

            # coefficients
            #c2 = 2                                                    # exponent
            #c1 = 1/24*(10**logg/10**(4.4))**(7/9)*(Teff/5777)**(23/9) # characteristic timescale
            #c0 = vpsd_avg[np.argmin(np.abs(freq_avg-1/c1))]           # amplitude

            # specifications
            name     = "Granulation"       # name
            type     = "Harvey"            # function type
            coef_val = [c0  , c1  , c2   ] # coefficient values
            coef_err = [0   , 0   , 0    ] # coefficient errors
            vary     = [True, True, False] # vary when fitted

            # check vpsd is computed
            if self.vpsd is not None:

                # check characteristic frequency is within frequency range
                if (1/c1 > min(freq)) & (1/c1 < max(freq)):

                    # save VPSD component
                    self.vpsd_components[name] = {"type": type, "coef_val": coef_val, "coef_err": coef_err, "vary": vary}
            
            else:

                # save VPSD component
                self.vpsd_components[name] = {"type": type, "coef_val": coef_val, "coef_err": coef_err, "vary": vary}

        # supergranulation
        if "Supergranulation" in components:

            # coefficients for the Sun
            c2_sun = 2                                    # exponent
            c1_sun = 1.33e1 / 24                          # characteristic timescale
            c0_sun = 2.79e4 * mps_to_kmps**2 * sec_to_day # amplitude

            # coefficients
            c2 = c2_sun                                                            # exponent
            c1 = c1_sun * (nu_max/nu_max_sun)**(-0.992)                            # characteristic timescale
            c0 = c0_sun * (nu_max/nu_max_sun)**(-0.609*2-0.992) * photo_to_velo**2 # amplitude

            # coefficients
            #c2 = 2                                                     # exponent
            #c1 = 10/24*(10**logg/10**(4.4))**(7/9)*(Teff/5777)**(23/9) # characteristic timescale
            #c0 = vpsd_avg[np.argmin(np.abs(freq_avg-1/c1))]            # amplitude

            # specifications
            name     = "Supergranulation"  # name
            type     = "Harvey"            # function type
            coef_val = [c0  , c1  , c2   ] # coefficient values
            coef_err = [0   , 0   , 0    ] # coefficient errors
            vary     = [True, True, False] # vary when fitted

            # check vpsd is computed
            if self.vpsd is not None:

                # check characteristic frequency is within frequency range
                if (1/c1 > min(freq)) & (1/c1 < max(freq)):

                    # save VPSD component
                    self.vpsd_components[name] = {"type": type, "coef_val": coef_val, "coef_err": coef_err, "vary": vary}
            
            else:

                # save VPSD component
                self.vpsd_components[name] = {"type": type, "coef_val": coef_val, "coef_err": coef_err, "vary": vary}

        # noise
        if "Noise" in components:

            # check vpsd is computed
            if self.vpsd is not None:

                # coefficients
                c0 = np.min(vpsd_avg) # amplitude

                # specifications
                name     = "Noise"    # name
                type     = "Constant" # function type
                coef_val = [c0  ]     # coefficient values
                coef_err = [0   ]     # coefficient errors
                vary     = [True]     # vary when fitted

                # save VPSD component
                self.vpsd_components[name] = {"type": type, "coef_val": coef_val, "coef_err": coef_err, "vary": vary}
        
        return None