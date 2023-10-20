import numpy             as     np
import pandas            as     pd
from   scipy.interpolate import interp1d
from   scipy.optimize    import curve_fit
from   tqdm              import tqdm
import warnings
warnings.filterwarnings('ignore')

class compute_vrad_pbp:

    def compute_vrad_pbp(self, mask_path:str=None, criteria:list=None, exclude_tellurics=True, scale=True, Niter=1, bins:"list[list]"=None) -> None:
        """Compute radial velocities (RVs) from spectral data.

        :param mask_path: path to line mask (must be a CSV file where the wavelength column is "wave"), defaults to None
        :type mask_path: str, optional
        :param criteria: criteria to apply (must be columns with prefix "crit_"), defaults to None
        :type criteria: list, optional
        :param exclude_tellurics: exclude telluric bands, defaults to True
        :type exclude_tellurics: bool, optional
        :param scale: scale master spectrum when fitting, defaults to True
        :type scale: bool, optional
        :param Niter: number of iterations per line, defaults to 1
        :type Niter: int, optional
        :param bins: temperature ranges with which to segment each line, defaults to None
        :type bins: list[list], optional
        :return: None
        :rtype: None
        """

        # read data
        wave_val                     =  self.spec["wave_val"]
        mast_wave_val, mast_flux_val = [self.spec_mast[key] for key in ["wave_val", "flux_val"]]
        temp                         =  self.aux_data["spec"]["temp"]
        if self.spec["path"] is None:
            flux_val_arr, flux_err_arr = [self.spec[key] for key in ["flux_val", "flux_err"]]

        # master spectrum gradient
        mast_grad_val = np.gradient(mast_flux_val)/np.gradient(mast_wave_val)
        
        # flux and gradient functions
        flux_func  = interp1d(mast_wave_val, mast_flux_val, kind="cubic")
        grad_func  = interp1d(mast_wave_val, mast_grad_val, kind="cubic")

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

        # temperature bins
        if bins is None:
            bins = [np.min(temp),np.max(temp)]

        # nr. of spectra, lines and temperature bins
        if self.spec["path"] is None:
            Nspec = len(flux_val_arr)
        else:
            Nspec = len(self.spec["files"])
        Nline = len(wc)
        Ntemp = len(bins)

        # lines which satisfy criteria
        idx_crit = np.ones_like(wc, dtype=bool)
        if criteria is not None:
            for i in range(len(criteria)):
                crit = np.array(mask["crit_"+criteria[i]])
                idx_crit *= crit

        # exclude tellurics
        idx_tell = np.ones_like(wc, dtype=bool)
        if exclude_tellurics:

            # read telluric bands
            tell      = self.aux_data["tell"]
            wave_l    = np.array(tell["wave_l"])
            wave_r    = np.array(tell["wave_r"])
            wave_band = np.array([wave_l,wave_r]).T

            # keep lines outside of bands
            idx_tell = np.sum([(wc > wave_band[i,0]) & (wc < wave_band[i,1]) for i in range(len(wave_band))], axis=0).astype(bool) == False

        # index of lines
        idx_line = np.zeros(Nline, dtype=object)
        for i in range(Nline):
            il          = np.searchsorted(wave_val, wl[i], side="right")
            ir          = np.searchsorted(wave_val, wr[i], side="right")
            idx_line[i] = np.arange(il, ir)

        # index of temperature bins
        idx_temp = np.zeros(Ntemp, dtype=object)
        for i in range(Ntemp):
            il          = temp>=bins[i][0]
            ir          = temp<=bins[i][1]
            idx_temp[i] = np.where(il & ir)[0]

        # index of temperature-segmented line parts
        idx = np.zeros((Nline,Ntemp), dtype=object)
        for i in range(Nline):
            for j in range(Ntemp):
                idx[i,j] = np.intersect1d(idx_line[i], idx_temp[j])

        # empty arrays for RV values and errors
        vrad_val_arr = np.zeros((Nspec,Nline,Ntemp))
        vrad_err_arr = np.zeros((Nspec,Nline,Ntemp))

        # loop spectra
        print("Analyzed spectra:")
        for i in tqdm(range(Nspec)):

            # read data from input
            if self.spec["path"] is None:
                flux_val = flux_val_arr[i]
                flux_err = flux_err_arr[i]
            
            # read data from path
            else:

                # read CSV files
                if self.spec["extension"] == "csv":

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
                    
                # read NPZ files
                if self.spec["extension"] == "npz":

                    # read NPZ file
                    file = np.load(self.spec["files"][i])
                    
                    # get flux and flux error if same wavelength grid
                    if self.spec["same_wave_grid"]:
                        flux_val   = file["flux_val"]
                        flux_err   = file["flux_err"]
                    
                    # interpolate flux and flux error on reference wavelength grid
                    else:
                        self.time["time_val"][i] = float(file["time_val"])
                        wave_val_i = file["wave_val"]
                        flux_val_i = file["flux_val"]
                        flux_err_i = file["flux_err"]
                        flux_val   = interp1d(wave_val_i, flux_val_i, kind="cubic", bounds_error=False)(wave_val)
                        flux_err   = interp1d(wave_val_i, flux_err_i, kind="cubic", bounds_error=False)(wave_val)

            # loop lines
            for j in range(Nline):

                # skip line if it does not satisfy criteria or overlaps with tellurics
                if idx_crit[j]*idx_tell[j] == False:

                    # assign NaN to RV and RV error
                    vrad_val_arr[i,j,:] = np.nan
                    vrad_err_arr[i,j,:] = np.nan
                
                # continue if good line
                else:

                    # loop temperature bins
                    for k in range(Ntemp):

                        # initial guess
                        vrad_val = 0

                        # wavelength, flux and flux uncertainty for line
                        wave_spec_line = wave_val[idx[j,k]]
                        flux_spec_line = flux_val[idx[j,k]]
                        ferr_spec_line = flux_err[idx[j,k]]

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
                            vrad_val_arr[i,j,k] = vrad_val
                            vrad_err_arr[i,j,k] = vrad_err
                        
                        # if unsuccessful, assign NaN to RV and RV error
                        except:

                            vrad_val_arr[i,j,k] = np.nan
                            vrad_err_arr[i,j,k] = np.nan

        # weighted average of valid lines
        vrad_val = np.zeros((Nspec,Ntemp))
        vrad_err = np.zeros((Nspec,Ntemp))
        for i in range(Nspec):
            for j in range(Ntemp):
                idx = ~np.isnan(vrad_val_arr[i,:,j]) & (np.abs(vrad_err_arr[i,:,j])>0)
                vrad_val[i,j] = np.average(vrad_val_arr[i,idx,j], weights=1/vrad_err_arr[i,idx,j]**2)
                vrad_err[i,j] = np.sqrt(1/np.sum(1/vrad_err_arr[i,idx,j]**2))
        
        # reduce dimensions if only one temperature bin
        if Ntemp == 1:
            vrad_val = vrad_val[:,0]
            vrad_err = vrad_err[:,0]
            vrad_val_arr = vrad_val_arr[:,:,0]
            vrad_err_arr = vrad_err_arr[:,:,0]

        # save RV data
        self.vrad = {"vrad_val": vrad_val, "vrad_err": vrad_err, "vrad_val_pbp": vrad_val_arr, "vrad_err_pbp": vrad_err_arr, "bins": bins, "method": "PBP", "mask": mask_name}

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