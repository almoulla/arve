import numpy             as     np
from   scipy.interpolate import interp1d
from   scipy.optimize    import curve_fit
from   tqdm              import tqdm

class compute_vrad_ccf:

    def compute_vrad_ccf(
        self,
        weight_name       : str         | None = None ,
        vrad_grid         : list[float] | None = None ,
        ccf_err_scale     : bool               = False,
        criteria          : list[str]   | None = None ,
        exclude_tellurics : bool               = True ,
        exclude_regions   : bool               = True ,
        ) -> None:
        """Compute radial velocities (RVs) from spectral data using the cross-correlation function (CCF) method.

        Parameters
        ----------
        weight_name : str | None, optional
            name of float column available in the mask to be used as weight (if None, all lines are weighted equally), by default None
        vrad_grid : list[float] | None, optional
            velocity grid, in the format [start,stop,step] and in units of km/s, on which to evaluate the CCF (if None, the grid is set to [-20,20,v_med] where v_med is the median wavelength step converted to a velocity), by default None
        ccf_err_scale : bool, optional
            scale the CCF errors if the velocity grid is under- or over-sampled, by default False
        criteria : list[str] | None, optional
            criteria to apply (must be boolean columns available in the mask), by default None
        exclude_tellurics : bool, optional
            exclude telluric bands, by default True
        exclude_regions : bool, optional
            exclude wavelength intervals, by default True

        Returns
        -------
        None
            None
        """

        # read data
        wave_val = self.spec["wave_val"]
        N_ord    = self.spec["N_ord"]

        # read constants
        c = self.arve.functions.constants["c"]

        # read mask
        mask = self.aux_data["mask"]

        # nr. of lines
        N_line = np.zeros(N_ord, dtype=int)
        for i in range(N_ord):
            N_line[i] = len(mask[i])

        # central wavelengths
        wave_c = np.zeros(N_ord, dtype=object)
        for i in range(N_ord):
            wave_c[i] = mask[i]["wave"].values

        # weights
        weight = np.zeros(N_ord, dtype=object)
        for i in range(N_ord):
            if weight_name is None:
                weight[i] = np.ones(N_line[i])
            else:
                weight[i] = mask[i][weight_name].values

        # exclude tellurics
        if exclude_tellurics:
            if criteria is None:
                criteria = ["crit_tell"]
            else:
                criteria.append("crit_tell")
        
        # exclude regions
        if exclude_regions:
            if criteria is None:
                criteria = ["crit_excl"]
            else:
                criteria.append("crit_excl")

        # lines which satisfy criteria
        if criteria is not None:
            for i in range(N_ord):
                idx = np.ones(N_line[i], dtype=bool)
                for j in range(len(criteria)):
                    crit = mask[i][criteria[j]].values
                    idx *= crit
                wave_c[i] = wave_c[i][idx]
                weight[i] = weight[i][idx]

        # RV shifts
        wave_val_flat = np.concatenate(wave_val)
        vrad_step_med = np.nanmedian(np.diff(wave_val_flat)/wave_val_flat[:-1])*c
        if vrad_grid is None:
            ccf_vrad_start = -20
            ccf_vrad_stop  =  20
            ccf_vrad_step  = vrad_step_med
            ccf_vrad       = np.arange(ccf_vrad_start, ccf_vrad_stop, ccf_vrad_step)
            ccf_err_corr   = 1
        else:
            ccf_vrad       = np.arange(vrad_grid[0], vrad_grid[1]+vrad_grid[2]/2, vrad_grid[2])
            ccf_err_corr   = (vrad_step_med/vrad_grid[2])**(1/2)

        # keep mask lines within spectrum overlap
        for i in range(N_ord):
            idx = (self.arve.functions.doppler_shift(wave=wave_c[i], v=np.min(ccf_vrad)) > np.min(wave_val[i])) & \
                  (self.arve.functions.doppler_shift(wave=wave_c[i], v=np.max(ccf_vrad)) < np.max(wave_val[i]))
            wave_c[i] = wave_c[i][idx]
            weight[i] = weight[i][idx]

        # normalize weights
        for i in range(N_ord):
            weight[i] = weight[i]/np.sum(weight[i])
        
        # nr. of spectra, RV shifts and lines
        if self.spec["path"] is None:
            N_spec = len(self.spec["flux_val"])
        else:
            N_spec = len(self.spec["files"])
        N_vrad = len(ccf_vrad)

        # empty arrays for indices
        i_l = np.zeros((N_ord,N_vrad), dtype=object)
        i_u = np.zeros((N_ord,N_vrad), dtype=object)

        # empty arrays for fractions
        r_l = np.zeros((N_ord,N_vrad), dtype=object)
        r_u = np.zeros((N_ord,N_vrad), dtype=object)

        # loop spectral orders
        for i in range(N_ord):

            # loop RV shifts
            for j in range(N_vrad):

                # shifted central wavelengths
                wc_shift = self.arve.functions.doppler_shift(wave=wave_c[i], v=ccf_vrad[j])

                # upper and lower indices
                i_u[i,j] = np.searchsorted(wave_val[i], wc_shift)
                i_l[i,j] = i_u[i,j] - 1

                # upper and lower fractions
                r_u[i,j] = (wc_shift-wave_val[i][i_l[i,j]])/(wave_val[i][i_u[i,j]]-wave_val[i][i_l[i,j]])
                r_l[i,j] = 1 - r_u[i,j]

        # empty arrays for CCF values and errors
        ccf_val = np.zeros((N_spec,N_ord+1,N_vrad))
        ccf_err = np.zeros((N_spec,N_ord+1,N_vrad))

        # empty arrays for CCF moments
        vrad_val_arr = np.zeros((N_spec,N_ord+1))*np.nan
        vrad_err_arr = np.zeros((N_spec,N_ord+1))*np.nan
        cont_val_arr = np.zeros((N_spec,N_ord+1))*np.nan
        cont_err_arr = np.zeros((N_spec,N_ord+1))*np.nan
        ampl_val_arr = np.zeros((N_spec,N_ord+1))*np.nan
        ampl_err_arr = np.zeros((N_spec,N_ord+1))*np.nan
        fwhm_val_arr = np.zeros((N_spec,N_ord+1))*np.nan
        fwhm_err_arr = np.zeros((N_spec,N_ord+1))*np.nan

        # loop spectra
        print("Analyzed spectra:")
        for i in tqdm(range(N_spec)):

            # read spectrum
            _, flux_val, flux_err = self.read_spec(i)

            # loop orders
            for j in range(N_ord+1):

                try:

                    # single order
                    if j < N_ord:

                        # loop RV shifts
                        for k in range(N_vrad):

                            # CCF value and error
                            ccf_val[i,j,k] = np.nansum(( flux_val[j][i_l[j,k]]*r_l[j,k]    + flux_val[j][i_u[j,k]]*r_u[j,k]    )*weight[j]   )
                            ccf_err[i,j,k] = np.nansum(((flux_err[j][i_l[j,k]]*r_l[j,k])**2+(flux_err[j][i_u[j,k]]*r_u[j,k])**2)*weight[j]**2)**(1/2)

                    # sum of all orders
                    else:

                        # order-summed CCF values and errors
                        ccf_val[i,j] = np.nansum(ccf_val[i]   , axis=0)
                        ccf_err[i,j] = np.nansum(ccf_err[i]**2, axis=0)**(1/2)
                    
                    # scale CCF errors
                    if ccf_err_scale:
                        ccf_err[i,j] *= ccf_err_corr

                    # initial guess on Gaussian parameters
                    i_min = np.argmin(ccf_val[i,j])
                    i_max = np.argmax(ccf_val[i,j])
                    C0 = ccf_val[i,j,i_max]
                    a0 = 1-ccf_val[i,j,i_min]/C0
                    b0 = ccf_vrad[i_min]
                    c0 = (b0-interp1d(ccf_val[i,j,:i_min], ccf_vrad[:i_min])((ccf_val[i,j,i_min]+ccf_val[i,j,i_max])/2))*2
                    if np.abs(b0) < 1e-3:
                        b0 = 1e-3
                    p0 = (C0, a0, b0, c0)

                    # CCF parameters with fitted Gaussian
                    param_val, param_cov = curve_fit(self.arve.functions.inverted_gaussian, ccf_vrad, ccf_val[i,j], sigma=ccf_err[i,j], p0=p0)
                    param_err            = np.sqrt(np.diag(param_cov))

                    # RV value and error
                    vrad_val_arr[i,j] = param_val[2]
                    vrad_err_arr[i,j] = 1/np.sqrt(np.sum(1/(ccf_err[i,j]*np.gradient(ccf_vrad)/np.gradient(ccf_val[i,j]))**2))

                    # continuum value and error
                    cont_val_arr[i,j] = param_val[0]
                    cont_err_arr[i,j] = param_err[0]

                    # amplitude value and error
                    ampl_val_arr[i,j] = param_val[1]
                    ampl_err_arr[i,j] = param_err[1]

                    # FWHM value and error
                    fwhm_val_arr[i,j] = param_val[3]
                    fwhm_err_arr[i,j] = param_err[3]
            
                except:
                    continue

        # emtpty arrays for weighted average of spectral orders
        vrad_val = np.zeros(N_spec)*np.nan
        vrad_err = np.zeros(N_spec)*np.nan
        cont_val = np.zeros(N_spec)*np.nan
        cont_err = np.zeros(N_spec)*np.nan
        ampl_val = np.zeros(N_spec)*np.nan
        ampl_err = np.zeros(N_spec)*np.nan
        fwhm_val = np.zeros(N_spec)*np.nan
        fwhm_err = np.zeros(N_spec)*np.nan

        # loop spectra
        for i in range(N_spec):

            try:

                # weighted average of RV
                vrad_val[i], vrad_err[i] = _weighted_average(vrad_val_arr[i,:-1], vrad_err_arr[i,:-1])
                
                # weighted average of continuum
                cont_val[i], cont_err[i] = _weighted_average(cont_val_arr[i,:-1], cont_err_arr[i,:-1])

                # weighted average of amplitude
                ampl_val[i], ampl_err[i] = _weighted_average(ampl_val_arr[i,:-1], ampl_err_arr[i,:-1])

                # weighted average of FWHM
                fwhm_val[i], fwhm_err[i] = _weighted_average(fwhm_val_arr[i,:-1], fwhm_err_arr[i,:-1])

            except:
                continue

        # bisector normalized flux points
        bis_flux_val = np.arange(0.00, 1.00, 0.01)
        N_bis        = len(bis_flux_val)

        # emtpty arrays for bisector and BIS value and error
        bis_vrad_val = np.zeros((N_spec,N_bis))*np.nan
        bis_vrad_err = np.zeros((N_spec,N_bis))*np.nan
        bis_val      = np.zeros( N_spec       )*np.nan
        bis_err      = np.zeros( N_spec       )*np.nan

        # loop spectra
        for i in range(N_spec):

            try:

                # CCF RV error
                ccf_err_vrad = (ccf_err[i,-1]*np.gradient(ccf_vrad)/np.gradient(ccf_val[i,-1]))**2

                # bisector value
                bis_vrad_val_l  = interp1d(ccf_val[i,-1][ccf_vrad<vrad_val[i]]/cont_val[i], ccf_vrad[ccf_vrad<vrad_val[i]], bounds_error=False)(bis_flux_val)
                bis_vrad_val_r  = interp1d(ccf_val[i,-1][ccf_vrad>vrad_val[i]]/cont_val[i], ccf_vrad[ccf_vrad>vrad_val[i]], bounds_error=False)(bis_flux_val)
                bis_vrad_val[i] = (bis_vrad_val_l+bis_vrad_val_r)/2

                # bisector error
                bis_vrad_err_l  = interp1d(ccf_vrad[ccf_vrad<vrad_val[i]], ccf_err_vrad[ccf_vrad<vrad_val[i]], bounds_error=False)(bis_vrad_val_l)
                bis_vrad_err_r  = interp1d(ccf_vrad[ccf_vrad>vrad_val[i]], ccf_err_vrad[ccf_vrad>vrad_val[i]], bounds_error=False)(bis_vrad_val_r)
                bis_vrad_err[i] = np.sqrt(bis_vrad_err_l**2 + bis_vrad_err_r**2)/2

                # bisector normalized flux
                bis_flux_norm = 1-(1-bis_flux_val)/ampl_val[i]

                # bisector index for bottom and top
                idx_bis_bot = (bis_flux_norm>0.1) & (bis_flux_norm<0.4)
                idx_bis_top = (bis_flux_norm>0.6) & (bis_flux_norm<0.9)

                # BIS bottom and top value and error
                bis_bot_val = np.average(bis_vrad_val[i][idx_bis_bot])
                bis_bot_err = np.average(bis_vrad_err[i][idx_bis_bot])/np.sqrt(np.sum(idx_bis_bot))
                bis_top_val = np.average(bis_vrad_val[i][idx_bis_top])
                bis_top_err = np.average(bis_vrad_err[i][idx_bis_top])/np.sqrt(np.sum(idx_bis_top))

                # BIS value and error
                bis_val[i] =  bis_top_val    - bis_bot_val
                bis_err[i] = (bis_top_err**2 + bis_bot_err**2)**(1/2)

            except:
                continue

        # save data
        self.vrad = {
            "vrad_val"    : vrad_val    ,
            "vrad_err"    : vrad_err    ,
            "vrad_val_ord": vrad_val_arr,
            "vrad_err_ord": vrad_err_arr,
            "method"      : "CCF"       ,
            }
        self.ccf  = {
            "ccf_vrad"    : ccf_vrad    ,
            "ccf_val"     : ccf_val     ,
            "ccf_err"     : ccf_err     ,
            "bis_vrad_val": bis_vrad_val,
            "bis_vrad_err": bis_vrad_err,
            "bis_flux_val": bis_flux_val,
            "cont_val"    : cont_val    ,
            "cont_err"    : cont_err    ,
            "ampl_val"    : ampl_val    ,
            "ampl_err"    : ampl_err    ,
            "fwhm_val"    : fwhm_val    ,
            "fwhm_err"    : fwhm_err    ,
            "bis_val"     : bis_val     ,
            "bis_err"     : bis_err     ,
            }

        return None

def _weighted_average(
    val_arr : np.ndarray,
    err_arr : np.ndarray
    ) -> tuple[float]:
    """Weighted average.

    Parameters
    ----------
    val_arr : np.ndarray
        value array
    err_arr : np.ndarray
        error array

    Returns
    -------
    tuple[float]
        weighted average value and error
    """

    idx = ~np.isnan(val_arr) & ~np.isnan(err_arr) & (err_arr>0)
    if np.sum(idx) > 0:
        val_avg = np.average(val_arr[idx], weights=1/err_arr[idx]**2)
        err_avg = np.sqrt(1/np.sum(1/err_arr[idx]**2))
    else:
        val_avg = np.nan
        err_avg = np.nan

    return val_avg, err_avg