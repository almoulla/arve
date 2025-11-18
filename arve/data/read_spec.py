from   astropy.io        import fits
import numpy             as     np
import pandas            as     pd
from   scipy.interpolate import interp1d

class read_spec:

    def read_spec(
        self,
        i : int,
        ) -> tuple[np.ndarray]:
        """Read spectrum.

        Parameters
        ----------
        i : int
            spectrum index

        Returns
        -------
        tuple[np.ndarray]
            wavelength values, flux values and flux errors of i:th spectrum
        """

        # read data from input
        if self.spec["files"] is None:
            wave_val = self.spec["wave_val"]
            flux_val = self.spec["flux_val"][i]
            flux_err = self.spec["flux_err"][i]
        
        # read data from files
        if self.spec["files"] is not None:

            # read file: FITS
            if self.spec["extension"] == "fits":

                # open file
                hdul = fits.open(self.spec["files"][i])

                # instrument: ESPRESSO
                if self.spec["instrument"] == "espresso":
                    if i == 0:
                        self.spec["medium"    ] = "air"
                        if hdul[0].header["HIERARCH ESO INS MODE"] == "SINGLEHR":
                            self.spec["resolution"] = 138000
                        if hdul[0].header["HIERARCH ESO INS MODE"] == "SINGLEUHR":
                            self.spec["resolution"] = 190000
                    self.time["time_val"][i] = hdul[0].header["HIERARCH ESO QC BJD"]
                    if self.spec["format"] == "s1d":
                        file = {"wave_val": hdul[1].data["wavelength_air"],
                                "flux_val": hdul[1].data["flux"          ],
                                "flux_err": hdul[1].data["error"         ]}
                    if self.spec["format"] == "s2d":
                        file = {"wave_val": hdul[5].data,
                                "flux_val": hdul[1].data,
                                "flux_err": hdul[2].data}
                
                # instrument: EXPRES
                if self.spec["instrument"] == "expres":
                    if i == 0:
                        self.spec["medium"    ] = "vac"
                        self.spec["resolution"] = 137500
                    self.time["time_val"][i] = hdul[1].header["BARYMJD"]
                    if self.spec["format"] == "s2d":
                        file = {"wave_val": hdul[1].data["bary_wavelength"]                      ,
                                "flux_val": hdul[1].data["spectrum"       ]*hdul[1].data["blaze"],
                                "flux_err": hdul[1].data["uncertainty"    ]*hdul[1].data["blaze"]}

                # instrument: HARPS
                if self.spec["instrument"] == "harps":
                    if i == 0:
                        self.spec["medium"    ] = "air"
                        self.spec["resolution"] = 115000
                    self.time["time_val"][i] = hdul[0].header["HIERARCH ESO QC BJD"]
                    if self.spec["format"] == "s1d":
                        file = {"wave_val": hdul[1].data["wavelength_air"],
                                "flux_val": hdul[1].data["flux"          ],
                                "flux_err": hdul[1].data["error"         ]}
                    if self.spec["format"] == "s2d":
                        file = {"wave_val": hdul[5].data,
                                "flux_val": hdul[1].data,
                                "flux_err": hdul[2].data}

                # instrument: HARPS-N
                if self.spec["instrument"] == "harps-n":
                    if i == 0:
                        self.spec["medium"    ] = "air"
                        self.spec["resolution"] = 115000
                    self.time["time_val"][i] = hdul[0].header["HIERARCH TNG QC BJD"]
                    if self.spec["format"] == "s1d":
                        file = {"wave_val": hdul[1].data["wavelength_air"],
                                "flux_val": hdul[1].data["flux"          ],
                                "flux_err": hdul[1].data["error"         ]}
                    if self.spec["format"] == "s2d":
                        file = {"wave_val": hdul[5].data,
                                "flux_val": hdul[1].data,
                                "flux_err": hdul[2].data}

                # instrument: NEID
                if self.spec["instrument"] == "neid":
                    if i == 0:
                        self.spec["medium"    ] = "vac"
                        self.spec["resolution"] = 110000
                    self.time["time_val"][i] = hdul[0].header["OBSJD"]
                    self.time["berv_val"][i] = hdul[0].header["SSBRV100"]*(-1)
                    if self.spec["format"] == "s2d":
                        file = {"wave_val": hdul[7].data       ,
                                "flux_val": hdul[1].data       ,
                                "flux_err": hdul[4].data**(1/2)}

                # instrument: NIRPS
                if self.spec["instrument"] == "nirps":
                    if i == 0:
                        self.spec["medium"    ] = "air"
                        self.spec["resolution"] = 80000
                    self.time["time_val"][i] = hdul[0].header["HIERARCH ESO QC BJD"]
                    if self.spec["format"] == "s1d":
                        file = {"wave_val": hdul[1].data["wavelength_air"],
                                "flux_val": hdul[1].data["flux"          ],
                                "flux_err": hdul[1].data["error"         ]}
                    if self.spec["format"] == "s2d":
                        file = {"wave_val": hdul[5].data,
                                "flux_val": hdul[1].data,
                                "flux_err": hdul[2].data}
                
                # instrument: SPIROU
                if self.spec["instrument"] == "spirou":
                    if i == 0:
                        self.spec["medium"    ] = "vac"
                        self.spec["resolution"] = 75000
                    self.time["time_val"][i] = hdul[0].header["MJDATE"]
                    self.time["berv_val"][i] = hdul[0].header["BERV"]
                    if self.spec["format"] == "s2d":
                        file = {"wave_val": hdul[1].data["wavelength"    ]*10     ,
                                "flux_val": hdul[1].data["flux"          ]        ,
                                "flux_err": hdul[1].data["weight"        ]**(-1/2)}
                
                # instrument: ESSP
                if self.spec["instrument"] == "essp":
                    if i == 0:
                        self.spec["medium"    ] = "vac"
                        self.spec["resolution"] = 100000
                    self.time["time_val"][i] = hdul[0].header["mjd_utc"]
                    if self.spec["format"] == "s2d":
                        file = {"wave_val": hdul["wavelength" ].data,
                                "flux_val": hdul["flux"       ].data,
                                "flux_err": hdul["uncertainty"].data}
                
                # close file
                hdul.close()
            
            # read file: NPZ
            if self.spec["extension"] == "npz":
                file = np.load(self.spec["files"][i])
                self.time["time_val"][i] = file["time_val"]
                try:
                    self.time["berv_val"][i] = file["berv_val"]
                except:
                    pass
            
            # read file: CSV
            if self.spec["extension"] == "csv":
                file = pd.read_csv(self.spec["files"][i])
            
            # reshape arrays
            wave_val = np.array(file["wave_val"]).astype("float64")
            flux_val = np.array(file["flux_val"]).astype("float64")
            flux_err = np.array(file["flux_err"]).astype("float64")
            if self.spec["format"] == "s1d":
                wave_val = wave_val.reshape(1,wave_val.shape[0])
                flux_val = flux_val.reshape(1,flux_val.shape[0])
                flux_err = flux_err.reshape(1,flux_err.shape[0])
            
            # replace points where wavelength is NaN or null with NaN
            idx_nan = np.isnan(wave_val) | (wave_val==0)
            wave_val[idx_nan] = np.nan
            flux_val[idx_nan] = np.nan
            flux_err[idx_nan] = np.nan
            
            # shift wavelengths
            vrad_sys = self.arve.star.stellar_parameters["vrad_sys"]
            berv_val = self.time["berv_val"]
            wave_val = self.arve.functions.doppler_shift(wave=wave_val, v=-vrad_sys)
            if self.spec["berv_corrected"] == False:
                wave_val = self.arve.functions.doppler_shift(wave=wave_val, v=-berv_val[i])

            # interpolate flux values and errors on reference wavelength grid
            if (self.spec["same_wave_grid"] == False) & (i > 0):
                kind           = self.spec["interpolation"]
                wave_val_inter = self.spec["wave_val"]
                flux_val_inter = np.zeros_like(wave_val_inter)
                flux_err_inter = np.zeros_like(wave_val_inter)
                for j in range(self.spec["N_ord"]):
                    idx_nan = np.isnan(wave_val[j]) | np.isnan(flux_val[j]) | np.isnan(flux_err[j])
                    idx_val = ~idx_nan
                    if np.sum(idx_val) > 0:
                        flux_val_inter[j,idx_val] = interp1d(wave_val[j,idx_val], flux_val[j,idx_val], kind=kind, bounds_error=False)(wave_val_inter[j,idx_val])
                        flux_err_inter[j,idx_val] = interp1d(wave_val[j,idx_val], flux_err[j,idx_val], kind=kind, bounds_error=False)(wave_val_inter[j,idx_val])
                    if np.sum(idx_nan) > 0:
                        flux_val_inter[j][idx_nan] = np.nan
                        flux_err_inter[j][idx_nan] = np.nan
                wave_val = wave_val_inter
                flux_val = flux_val_inter
                flux_err = flux_err_inter
        
        return wave_val, flux_val, flux_err