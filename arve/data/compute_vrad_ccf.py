import numpy             as     np
import pandas            as     pd
from   scipy.interpolate import interp1d
from   scipy.optimize    import curve_fit
from   tqdm              import tqdm

class compute_vrad_ccf:

    def compute_vrad_ccf(self, mask_path:str=None, weight:str=None, criteria:list=None, exclude_tellurics=True, vgrid:list=[-20,20,0.25]) -> None:
        """Compute radial velocities (RVs) from spectral data.

        :param mask_path: path to line mask (must be a CSV file where the wavelength column is "wave"), defaults to None
        :type mask_path: str, optional
        :param weight: column name of weight, defaults to None
        :type weight: str, optional
        :param criteria: criteria to apply (must be columns with prefix "crit_"), defaults to None
        :type criteria: list, optional
        :param exclude_tellurics: exclude telluric bands, defaults to True
        :type exclude_tellurics: bool, optional
        :param vgrid: velocity grid, in the format [start,stop,step] and in units of km/s, on which to evaluate the CCF, defaults to [-20,20,0.25]
        :type vgrid: list, optional
        :return: None
        :rtype: None
        """

        # read data
        wave_val = self.spec["wave_val"]
        if self.spec["path"] is None:
            flux_val_arr, flux_err_arr = [self.spec[key] for key in ["flux_val", "flux_err"]]

        # mask from VALD
        if mask_path is None:

            # read mask and mask name
            mask      = self.aux_data["mask"]
            mask_name = self.aux_data["name"]
        
        # provided mask
        else:

            # read mask and mask name
            mask      = pd.read_csv(mask_path)
            mask_name = mask_path.split("/")[-1]

        # central wavelengths
        wc = np.array(mask["wave"])

        # weights
        if weight is None:
            w = np.ones_like(wc)
        else:
            w = np.array(mask[weight])

        # keep mask lines which satisfy criteria
        if criteria is not None:
            idx = np.ones_like(wc, dtype=bool)
            for i in range(len(criteria)):
                crit = np.array(mask["crit_"+criteria[i]])
                idx *= crit
            wc  = wc[idx]
            w   = w [idx]

        # exclude tellurics
        if exclude_tellurics:

            # read telluric bands
            tell      = self.aux_data["tell"]
            wave_l    = np.array(tell["wave_l"])
            wave_r    = np.array(tell["wave_r"])
            wave_band = np.array([wave_l,wave_r]).T

            # keep lines outside of bands
            idx = np.sum([(wc > wave_band[i,0]) & (wc < wave_band[i,1]) for i in range(len(wave_band))], axis=0).astype(bool) == False
            wc  = wc[idx]
            w   = w [idx]

        # RV shifts
        vrads = np.arange(vgrid[0],vgrid[1]+vgrid[2]/2,vgrid[2])

        # keep mask lines within spectrum overlap
        idx = (self.arve.functions.doppler_shift(wave=wc, v=min(vrads)) > min(wave_val)) & \
              (self.arve.functions.doppler_shift(wave=wc, v=max(vrads)) < max(wave_val))
        wc  = wc[idx]
        w   = w [idx]

        # normalize weights
        w = w/np.sum(w)
        
        # nr. of spectra
        if self.spec["path"] is None:
            Nspec = len(flux_val_arr)
        else:
            Nspec = len(self.spec["files"])
        
        # nr. of RV shifts
        Nvrad = len(vrads)
        
        # nr. of lines
        Nline = len(wc)

        # empty arrays for indices
        ir = np.zeros((Nvrad,Nline), dtype=int)
        il = np.zeros((Nvrad,Nline), dtype=int)

        # empty arrays for fractions
        fr = np.zeros((Nvrad,Nline))
        fl = np.zeros((Nvrad,Nline))

        # loop RV shifts
        for i in range(Nvrad):

            # shifted central wavelengths
            wc_shift = self.arve.functions.doppler_shift(wave=wc, v=vrads[i])

            # right and left indices
            ir[i] = np.searchsorted(wave_val, wc_shift, side="right")
            il[i] = ir[i] - 1

            # right and left fractions
            fr[i] = (wc_shift-wave_val[il[i]])/(wave_val[ir[i]]-wave_val[il[i]])
            fl[i] = 1 - fr[i]

        # empty arrays for RV values and errors
        vrad_val = np.zeros(Nspec)
        vrad_err = np.zeros(Nspec)

        # empty arrays for FWHM values and errors
        fwhm_val = np.zeros(Nspec)

        # loop spectra
        print("Analyzed spectra:")
        for i in tqdm(range(Nspec)):

            # read data from input
            if self.spec["path"] is None:
                flux_val = flux_val_arr[i]
                flux_err = flux_err_arr[i]
            
            # read data from path
            else:

                # read CSV file
                df = pd.read_csv(self.spec["files"][i])
                
                # get flux and flux error if same wavelength grid
                if self.spec["same_wave_grid"]:
                    flux_val   = df["flux_val"].to_numpy()
                    flux_err   = df["flux_err"].to_numpy()
                
                # interpolate flux and flux error on reference wavelength grid
                else:
                    wave_val_i = df["wave_val"].to_numpy()
                    flux_val_i = df["flux_val"].to_numpy()
                    flux_err_i = df["flux_err"].to_numpy()
                    flux_val   = interp1d(wave_val_i, flux_val_i, kind="cubic", bounds_error=False)(wave_val)
                    flux_err   = interp1d(wave_val_i, flux_err_i, kind="cubic", bounds_error=False)(wave_val)

            # empty arrays for CCF values and errors
            ccf_val = np.zeros(Nvrad)
            ccf_err = np.zeros(Nvrad)

            # loop RV shifts
            for j in range(Nvrad):

                # CCF value and error
                ccf_val[j] = np.nansum((flux_val[il[j]]   *fl[j]+flux_val[ir[j]]   *fr[j])*w   )
                ccf_err[j] = np.nansum((flux_err[il[j]]**2*fl[j]+flux_err[ir[j]]**2*fr[j])*w**2)**(1/2)

            # initial guess on Gaussian parameters
            i_min = np.argmin(ccf_val)
            i_max = np.argmax(ccf_val)
            C0 = ccf_val[i_max]
            a0 = 1-ccf_val[i_min]/C0
            b0 = vrads[i_min]
            c0 = (b0-interp1d(ccf_val[:i_min], vrads[:i_min], kind="cubic")((ccf_val[i_min]+ccf_val[i_max])/2))*2
            p0 = (C0, a0, b0, c0)

            # CCF parameters with fitted Gaussian
            param, _ = curve_fit(self.arve.functions.inverted_gaussian, vrads, ccf_val, sigma=ccf_err, p0=p0)

            # RV value and error
            vrad_val[i] = param[2]
            vrad_err[i] = 1/np.sqrt(np.sum(1/np.abs(ccf_err*np.gradient(vrads)/np.gradient(ccf_val))**2))

            # FWHM value and error
            fwhm_val[i] = param[3]

        # save RV data
        self.vrad = {"vrad_val": vrad_val, "vrad_err": vrad_err, "method": "CCF", "mask": mask_name}
        self.ccf  = {"fwhm_val": fwhm_val}

        return None