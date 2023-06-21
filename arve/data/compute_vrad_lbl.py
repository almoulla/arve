import numpy             as     np
import pandas            as     pd
from   scipy.interpolate import interp1d
from   scipy.optimize    import curve_fit
from   tqdm              import tqdm
import warnings
warnings.filterwarnings('ignore')

class compute_vrad_lbl:

    def compute_vrad_lbl(self, mask_path:str=None, weight:str=None, criteria:list=None, exclude_tellurics=True, scale=True, Niter=1) -> None:
        """Compute radial velocities (RVs) from spectral data.

        :param mask_path: path to line mask (must be a CSV file where the wavelength column is "wave"), defaults to None
        :type mask_path: str, optional
        :param weight: column name of weight, defaults to None
        :type weight: str, optional
        :param criteria: criteria to apply (must be columns with prefix "crit_"), defaults to None
        :type criteria: list, optional
        :param exclude_tellurics: exclude telluric bands, defaults to True
        :type exclude_tellurics: bool, optional
        :param scale: scale master spectrum when fitting, defaults to True
        :type scale: bool, optional
        :param Niter: number of iterations per line, defaults to 1
        :type Niter: int, optional
        :return: None
        :rtype: None
        """

        # read data
        wave_val                     =  self.spec["wave_val"]
        mast_wave_val, mast_flux_val = [self.spec_mast[key] for key in ["wave_val", "flux_val"]]
        if self.spec["path"] is None:
            flux_val_arr, flux_err_arr = [self.spec[key] for key in ["flux_val", "flux_err"]]

        # read constants
        c = self.arve.functions.constants["c"]

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

        # central, left and right wavelengths
        wc = np.array(mask["wave"  ])
        wl = np.array(mask["wave_l"])
        wr = np.array(mask["wave_r"])

        # keep mask lines which satisfy criteria
        if criteria is not None:
            idx = np.ones_like(wc, dtype=bool)
            for i in range(len(criteria)):
                crit = np.array(mask["crit_"+criteria[i]])
                idx *= crit
            wc  = wc[idx]
            wl  = wl[idx]
            wr  = wr[idx]

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
            wl  = wl[idx]
            wr  = wr[idx]
        
        # nr. of spectra
        if self.spec["path"] is None:
            Nspec = len(flux_val_arr)
        else:
            Nspec = len(self.spec["files"])

        # nr. of lines
        Nline = len(wc)

        # empty arrays for RV values and errors
        vrad_val_arr = np.zeros((Nspec,Nline))
        vrad_err_arr = np.zeros((Nspec,Nline))

        # master spectrum gradient
        mast_grad_val = np.gradient(mast_flux_val)/np.gradient(mast_wave_val)
        
        # flux and gradient functions
        flux_func  = interp1d(mast_wave_val, mast_flux_val, kind="cubic")
        grad_func  = interp1d(mast_wave_val, mast_grad_val, kind="cubic")

        # empty array for line indices
        ilr = np.zeros(Nline, dtype=object)

        # loop lines
        for i in range(Nline):
        
            # line index range in wavelength grid
            il     = np.searchsorted(wave_val, wl[i], side="right")
            ir     = np.searchsorted(wave_val, wr[i], side="right")
            ilr[i] = np.arange(il, ir)

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

            # loop lines
            for j in range(Nline):

                # initial guess
                vrad_val = 0

                # wavelength, flux and flux uncertainty for line
                wave_spec_line = wave_val[ilr[j]]
                flux_spec_line = flux_val[ilr[j]]
                ferr_spec_line = flux_err[ilr[j]]

                # try to template match
                try:

                    # loop iterations
                    for _ in range(Niter):

                        # single spectrum shift
                        wave_spec_shift = wave_spec_line*(1-vrad_val/c)
                        flux_spec_shift = flux_spec_line
                        ferr_spec_shift = ferr_spec_line
                        
                        # master spectrum interpolate
                        wave_mast_inter = wave_spec_shift
                        flux_mast_inter = flux_func(wave_mast_inter)
                        grad_mast_inter = grad_func(wave_mast_inter)

                        # select function of shifted spectrum and initial guess on parameters
                        if scale:
                            func_spec_shift = _spec_shift_scale
                            p0              = (0,0,1)
                        else:
                            func_spec_shift = _spec_shift
                            p0              = (0,)

                        # solve for parameters
                        param, _ = curve_fit(func_spec_shift, (wave_mast_inter, flux_mast_inter, grad_mast_inter), flux_spec_shift, sigma=ferr_spec_shift, p0=p0)

                        # compute RV and RV error
                        vrad_val += param[0]*c
                        vrad_err  = 1/np.sqrt(np.sum(1/(ferr_spec_shift/(grad_mast_inter*wave_mast_inter/c))**2))

                    # save
                    vrad_val_arr[i,j] = vrad_val
                    vrad_err_arr[i,j] = vrad_err
                
                # if unsuccessful, assign NaN to RV and RV error
                except:

                    vrad_val_arr[i,j] = np.nan
                    vrad_err_arr[i,j] = np.nan

        # weighted average of valid lines
        vrad_val = np.zeros(Nspec)
        vrad_err = np.zeros(Nspec)
        for i in range(Nspec):
            idx = ~np.isnan(vrad_val_arr[i,:]) & (np.abs(vrad_err_arr[i,:])>0)
            vrad_val[i] = np.average(vrad_val_arr[i,idx], weights=1/vrad_err_arr[i,idx]**2)
            vrad_err[i] = np.sqrt(1/np.sum(1/vrad_err_arr[i,idx]**2))

        # save RV data
        self.vrad = {"vrad_val": vrad_val, "vrad_err": vrad_err, "vrad_val_lbl": vrad_val_arr, "vrad_err_lbl": vrad_err_arr, "method": "LBL", "mask": mask_name}

        return None

# spectrum shifted
def _spec_shift(X, *param):

    wave, flux, grad = X
    z = param
    S = flux - grad*wave*z

    return S

# spectrum shfited and scaled
def _spec_shift_scale(X, *param):

    wave, flux, grad = X
    z, c0, c1 = param
    S = c0 + c1*(flux - grad*wave*z)

    return S