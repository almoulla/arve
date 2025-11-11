import numpy as np

from typing import Literal

class add_vpsd_components:

    def add_vpsd_components(
        self,
        components : list[Literal["oscillations", "granulation", "supergranulation", "noise"]] = ["oscillations", "granulation", "supergranulation", "noise"],
        ) -> None:
        """Add velocity power spectral density (VPSD) components.

        Parameters
        ----------
        components : list[Literal["oscillations", "granulation", "supergranulation", "noise"]], optional
            VPSD components, by default ["oscillations", "granulation", "supergranulation", "noise"]

        Returns
        -------
        None
            None
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
        microHz_to_Hz     = 1e-6                                        # micro-Hertz        to Hertz
        milliHz_to_Hz     = 1e-3                                        # milli-Hertz        to Hertz
        mps_to_kmps       = 1e-3                                        # m/s                to km/s
        sec_to_day        = 1/(60*60*24)                                # seconds            to days
        hr_to_day         = 1/(24)                                      # hours              to days
        std_to_hwhm       = np.sqrt(2*np.log(2))                        # standard deviation to half-width at half-maxmimum (HWHM)
        
        # conversions
        photo_to_velo_osc = (Teff/Teff_sun)**(1.8)                      # photometry         to velocimetry scaling ratio for oscillations
        photo_to_velo_gra = (Teff/Teff_sun)**(64/9) * (g/g_sun)**(-4/9) # photometry         to velocimetry scaling ratio for (super)granulation

        # frequency of maximum power (central oscillation mode)
        freq_max_sun = 2.73e+2 # 1/(3.66e-3)
        freq_max     = freq_max_sun * (Teff/Teff_sun)**(-1/2) * (g/g_sun)

        # dictionary with VPSD components
        if self.vpsd_components is None:
            self.vpsd_components = {}

        # oscillations
        if "oscillations" in components:

            # coefficients for the Sun
            c2_sun = freq_max_sun # central frequency, 1/(3.66e-3)
            c1_sun = 2.50e+1      # envelope HWHM    , 2.89e-1 * milliHz_to_Hz / sec_to_day
            c0_sun = 1.45e-8      # amplitude        , 1.25e+3 * mps_to_kmps**2 * sec_to_day

            # coefficients
            c2 = freq_max                                                       # central frequency,
            c1 = c1_sun * (freq_max/freq_max_sun)**( 0.880)                     # envelope HWHM    ,
            c0 = c0_sun * (freq_max/freq_max_sun)**(-2.305) * photo_to_velo_osc # amplitude        , -0.609*2/0.855-0.880 = -2.305

            # specifications
            comp_name = "oscillations"        # component name
            func_type = "lorentz"             # function type
            coef_name = ["amplitude"        , # coefficient names
                         "envelope hwhm"    ,
                         "central frequency"]
            coef_val  = [c0  , c1  , c2  ]    # coefficient values
            coef_err  = [0   , 0   , 0   ]    # coefficient errors
            coef_vary = [True, True, True]    # vary coefficients when fitted

            # check VPSD is computed
            if self.vpsd is not None:

                # check central frequency is within frequency range
                if (c2 > min(freq)) & (c2 < max(freq)):

                    # save VPSD component
                    self.vpsd_components[comp_name] = {"func_type": func_type, "coef_name": coef_name, "coef_val": coef_val, "coef_err": coef_err, "coef_vary": coef_vary}
            
            else:

                # save VPSD component
                self.vpsd_components[comp_name] = {"func_type": func_type, "coef_name": coef_name, "coef_val": coef_val, "coef_err": coef_err, "coef_vary": coef_vary}

        # granulation
        if "granulation" in components:

            # coefficients for the Sun
            c2_sun = 2       # decay rate,
            c1_sun = 3.80e-2 # timescale , 9.11e-1 * hr_to_day
            c0_sun = 4.76e-9 # amplitude , 4.11e+2 * mps_to_kmps**2 * sec_to_day

            # coefficients
            c2 = c2_sun                                                         # decay rate,
            c1 = c1_sun * (freq_max/freq_max_sun)**(-0.970)                     # timescale ,
            c0 = c0_sun * (freq_max/freq_max_sun)**(-2.188) * photo_to_velo_gra # amplitude , -0.609*2-0.970 = -2.188

            # specifications
            comp_name = "granulation"         # component name
            func_type = "harvey"              # function type
            coef_name = ["amplitude" ,        # coefficient names
                         "timescale" ,
                         "decay rate"]
            coef_val  = [c0  , c1  , c2   ]   # coefficient values
            coef_err  = [0   , 0   , 0    ]   # coefficient errors
            coef_vary = [True, True, False]   # vary coefficients when fitted

            # check VPSD is computed
            if self.vpsd is not None:

                # check characteristic frequency is within frequency range
                if (1/c1 > min(freq)) & (1/c1 < max(freq)):

                    # save VPSD component
                    self.vpsd_components[comp_name] = {"func_type": func_type, "coef_name": coef_name, "coef_val": coef_val, "coef_err": coef_err, "coef_vary": coef_vary}
            
            else:

                # save VPSD component
                self.vpsd_components[comp_name] = {"func_type": func_type, "coef_name": coef_name, "coef_val": coef_val, "coef_err": coef_err, "coef_vary": coef_vary}

        # supergranulation
        if "supergranulation" in components:

            # coefficients for the Sun
            c2_sun = 2       # decay rate,
            c1_sun = 5.54e-1 # timescale , 1.33e+1 * hr_to_day
            c0_sun = 3.23e-7 # amplitude , 2.79e+4 * mps_to_kmps**2 * sec_to_day

            # coefficients
            c2 = c2_sun                                                         # decay rate,
            c1 = c1_sun * (freq_max/freq_max_sun)**(-0.992)                     # timescale ,
            c0 = c0_sun * (freq_max/freq_max_sun)**(-2.210) * photo_to_velo_gra # amplitude , -0.609*2-0.992 = -2.210

            # specifications
            comp_name = "supergranulation"    # name
            func_type = "harvey"              # function type
            coef_name = ["amplitude" ,        # coefficient names
                         "timescale" ,
                         "decay rate"]
            coef_val  = [c0  , c1  , c2   ]   # coefficient values
            coef_err  = [0   , 0   , 0    ]   # coefficient errors
            coef_vary = [True, True, False]   # vary coefficients when fitted

            # check VPSD is computed
            if self.vpsd is not None:

                # check characteristic frequency is within frequency range
                if (1/c1 > min(freq)) & (1/c1 < max(freq)):

                    # save VPSD component
                    self.vpsd_components[comp_name] = {"func_type": func_type, "coef_name": coef_name, "coef_val": coef_val, "coef_err": coef_err, "coef_vary": coef_vary}
            
            else:

                # save VPSD component
                self.vpsd_components[comp_name] = {"func_type": func_type, "coef_name": coef_name, "coef_val": coef_val, "coef_err": coef_err, "coef_vary": coef_vary}

        # noise
        if "noise" in components:

            # check VPSD is computed
            if self.vpsd is not None:

                # coefficients
                c0 = np.min(vpsd_avg) # amplitude

                # specifications
                comp_name = "noise"       # component name
                func_type = "constant"    # function type
                coef_name = ["amplitude"] # coefficient names
                coef_val  = [c0  ]        # coefficient values
                coef_err  = [0   ]        # coefficient errors
                coef_vary = [True]        # vary coefficients when fitted

                # save VPSD component
                self.vpsd_components[comp_name] = {"func_type": func_type, "coef_name": coef_name, "coef_val": coef_val, "coef_err": coef_err, "coef_vary": coef_vary}
        
        return None