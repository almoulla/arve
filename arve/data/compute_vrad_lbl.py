from   astropy.stats     import sigma_clip
import numpy             as     np
from   scipy.interpolate import interp1d
from   scipy.optimize    import curve_fit
from   tqdm              import tqdm
import warnings
warnings.filterwarnings("ignore")

class compute_vrad_lbl:

    def compute_vrad_lbl(
        self,
        scale             : bool                     = True,
        bins              : list[list[float]] | None = None,
        N_iter            : int                      = 1   ,
        vrad_err_lim      : float                    = 1e-3,
        criteria          : list[str]         | None = None,
        exclude_tellurics : bool                     = True,
        exclude_regions   : bool                     = True,
        ) -> None:
        """Compute radial velocities (RVs) from spectral data using the line-by-line (LBL) method.

        Parameters
        ----------
        scale : bool, optional
            scale reference spectrum when fitting, by default True
        bins : list[list[float]] | None, optional
            formation temperature bins with which to segment each line (if None, only 1 bin with the entire line-forming temperature range is considered), by default None
        N_iter : int, optional
            number of iterations per line segment (the inital RV guess is updated from the previous iteration), by default 1
        vrad_err_lim : float, optional
            RV error lower limit (used to prevent diverging RVs due to null errors returned by the numerical solver for some lines), by default 1e-3
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
        wave_val                   =  self.spec["wave_val"]
        ref_wave_val, ref_flux_val = [self.spec_reference[key] for key in ["wave_val", "flux_val"]]
        temp                       =  self.aux_data["spec"]["temp"]
        N_spec                     =  self.spec["N_spec"]
        N_ord                      =  self.spec["N_ord"]

        # reference spectrum gradient
        ref_grad_val = np.zeros_like(ref_wave_val)
        for i in range(N_ord):
            ref_grad_val[i] = np.gradient(ref_flux_val[i])/np.gradient(ref_wave_val[i])
        
        # reference flux and gradient functions
        ref_flux_val_func = np.zeros(N_ord, dtype=object)
        ref_grad_val_func = np.zeros(N_ord, dtype=object)
        for i in range(N_ord):
            ref_flux_val_func[i] = interp1d(ref_wave_val[i], ref_flux_val[i], kind="cubic")
            ref_grad_val_func[i] = interp1d(ref_wave_val[i], ref_grad_val[i], kind="cubic")

        # read constants
        c = self.arve.functions.constants["c"]

        # read mask
        mask = self.aux_data["mask"]

        # temperature bins
        if bins is None:
            bins = [[np.nanmin(temp),np.nanmax(temp)]]

        # nr. of lines and temperature bins
        N_line = np.zeros(N_ord, dtype=int)
        for i in range(N_ord):
            N_line[i] = len(mask[i])
        N_bin = len(bins)

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
        idx_crit = np.zeros(N_ord, dtype=object)
        for i in range(N_ord):
            idx_crit[i] = np.ones(N_line[i], dtype=bool)
            if criteria is not None:
                for j in range(len(criteria)):
                    crit = np.array(mask[i][criteria[j]])
                    idx_crit[i] *= crit

        # lower and upper index of lines
        idx_l = np.zeros(N_ord, dtype=object)
        idx_u = np.zeros(N_ord, dtype=object)
        for i in range(N_ord):
            idx_l[i] = np.array(mask[i]["idx_l"])
            idx_u[i] = np.array(mask[i]["idx_u"])

        # index of lines
        idx_line = np.zeros(N_ord, dtype=object)
        for i in range(N_ord):
            idx_line[i] = np.zeros(N_line[i], dtype=object)
            for j in range(N_line[i]):
                idx_line[i][j] = np.arange(idx_l[i][j], idx_u[i][j]+1)

        # index of temperature bins
        idx_temp = np.zeros(N_ord, dtype=object)
        for i in range(N_ord):
            idx_temp[i] = np.zeros(N_bin, dtype=object)
            for j in range(N_bin):
                idx_temp[i][j] = np.where((temp[i]>=bins[j][0]) & (temp[i]<=bins[j][1]))[0]

        # index of temperature-segmented line parts
        idx = np.zeros(N_ord, dtype=object)
        for i in range(N_ord):
            idx[i] = np.zeros((N_line[i],N_bin), dtype=object)
            for j in range(N_line[i]):
                for k in range(N_bin):
                    idx[i][j,k] = np.intersect1d(idx_line[i][j], idx_temp[i][k])

        # empty arrays for RV values and errors
        vrad_val_lbl = np.zeros(N_ord, dtype=object)
        vrad_err_lbl = np.zeros(N_ord, dtype=object)
        for i in range(N_ord):
            vrad_val_lbl[i] = np.zeros((N_spec,N_line[i],N_bin))*np.nan
            vrad_err_lbl[i] = np.zeros((N_spec,N_line[i],N_bin))*np.nan
        
        # mask array for RV values and errors
        vrad_mask_lbl = np.zeros(N_ord, dtype=object)
        for i in range(N_ord):
            vrad_mask_lbl[i] = np.ones((N_spec,N_line[i],N_bin), dtype=bool)

        # loop spectra
        print("Analyzed spectra:")
        for i in tqdm(range(N_spec)):

            # read spectrum
            _, flux_val, flux_err = self.read_spec(i)

            # loop orders
            for j in range(N_ord):

                # loop lines
                for k in range(N_line[j]):

                    # check if line satisfies criteria
                    if idx_crit[j][k]:

                        # loop temperature bins
                        for l in range(N_bin):

                            # initial guess
                            vrad_val = 0

                            # wavelength, flux and flux error for line part
                            obs_wave_val = wave_val[j,idx[j][k,l]]
                            obs_flux_val = flux_val[j,idx[j][k,l]]
                            obs_flux_err = flux_err[j,idx[j][k,l]]

                            # try to template match
                            try:

                                # loop iterations
                                for _ in range(N_iter):

                                    # observed spectrum shifted
                                    obs_wave_val_shift = obs_wave_val*(1-vrad_val/c)
                                    obs_flux_val_shift = obs_flux_val
                                    obs_flux_err_shift = obs_flux_err
                                    
                                    # reference spectrum interpolated
                                    ref_wave_val_inter = obs_wave_val_shift
                                    ref_flux_val_inter = ref_flux_val_func[j](ref_wave_val_inter)
                                    ref_grad_val_inter = ref_grad_val_func[j](ref_wave_val_inter)

                                    # select function of shifted spectrum and initial guess on parameters
                                    if scale:
                                        func_spec_shift = _spec_shift_scale
                                        p0              = (0,1)
                                    else:
                                        func_spec_shift = _spec_shift
                                        p0              = (0,)

                                    # solve for parameters
                                    param, _ = curve_fit(func_spec_shift, (ref_wave_val_inter, ref_flux_val_inter, ref_grad_val_inter), obs_flux_val_shift, sigma=obs_flux_err_shift, p0=p0)

                                    # compute RV and RV error
                                    vrad_val += param[0]*c
                                    vrad_err  = 1/np.sqrt(np.sum(1/(obs_flux_err_shift/(ref_grad_val_inter*ref_wave_val_inter/c))**2))
                                    if scale:
                                        vrad_err /= np.sqrt(param[1])

                                # save
                                vrad_val_lbl[j][i,k,l] = vrad_val
                                vrad_err_lbl[j][i,k,l] = vrad_err
                            
                            # if unsuccessful, continue
                            except:

                                continue

        # replace infinite values with NaNs
        for i in range(N_ord):
            vrad_val_lbl[i][np.abs(vrad_val_lbl[i]) == np.inf] = np.nan
            vrad_err_lbl[i][np.abs(vrad_err_lbl[i]) == np.inf] = np.nan

        # update mask array
        for i in range(N_ord):
            vrad_mask_lbl[i][np.isnan(vrad_val_lbl[i]) | np.isnan(vrad_err_lbl[i])] = False

        # outlier identification
        for i in range(N_spec):
            *_, clip_val_min, clip_val_max = sigma_clip(np.vstack([vrad_val_lbl[j][i,:,:] for j in range(N_ord)]), maxiters=None, return_bounds=True)
            *_, clip_err_min, clip_err_max = sigma_clip(np.vstack([vrad_err_lbl[j][i,:,:] for j in range(N_ord)]), maxiters=None, return_bounds=True)
            for j in range(N_ord):
                for k in range(N_line[j]):
                    for l in range(N_bin):
                        if (vrad_val_lbl[j][i,k,l] < clip_val_min) | (vrad_val_lbl[j][i,k,l] > clip_val_max):
                            vrad_mask_lbl[j][i,k,l] = False
                        if (vrad_err_lbl[j][i,k,l] < clip_err_min) | (vrad_err_lbl[j][i,k,l] > clip_err_max):
                            vrad_mask_lbl[j][i,k,l] = False

        # weighted average of all valid lines per order
        vrad_val_ord = np.zeros((N_spec,N_ord,N_bin))*np.nan
        vrad_err_ord = np.zeros((N_spec,N_ord,N_bin))*np.nan
        for i in range(N_spec):
            for j in range(N_ord):
                for k in range(N_bin):
                    idx = vrad_mask_lbl[j][i,:,k] & (vrad_err_lbl[j][i,:,k]>vrad_err_lim)
                    if np.sum(idx) > 0:
                        vrad_val_ord[i,j,k] = np.average(vrad_val_lbl[j][i,idx,k], weights=1/vrad_err_lbl[j][i,idx,k]**2)
                        vrad_err_ord[i,j,k] = np.sqrt(1/np.sum(1/vrad_err_lbl[j][i,idx,k]**2))

        # weighted average of all orders per temperature bin
        vrad_val_bin = np.zeros((N_spec,N_bin))*np.nan
        vrad_err_bin = np.zeros((N_spec,N_bin))*np.nan
        for i in range(N_spec):
            for j in range(N_bin):
                idx = ~np.isnan(vrad_val_ord[i,:,j]) & (vrad_err_ord[i,:,j]>0)
                if np.sum(idx) > 0:
                    vrad_val_bin[i,j] = np.average(vrad_val_ord[i,idx,j], weights=1/vrad_err_ord[i,idx,j]**2)
                    vrad_err_bin[i,j] = np.sqrt(1/np.sum(1/vrad_err_ord[i,idx,j]**2))

        # default RV time series (taken to be the first temperature bin)
        vrad_val = vrad_val_bin[:,0]
        vrad_err = vrad_err_bin[:,0]

        # save RV data
        self.vrad = {
            "vrad_val"     : vrad_val     ,
            "vrad_err"     : vrad_err     ,
            "vrad_val_bin" : vrad_val_bin ,
            "vrad_err_bin" : vrad_err_bin ,
            "vrad_val_ord" : vrad_val_ord ,
            "vrad_err_ord" : vrad_err_ord ,
            "vrad_val_lbl" : vrad_val_lbl ,
            "vrad_err_lbl" : vrad_err_lbl ,
            "vrad_mask_lbl": vrad_mask_lbl,
            "bins"         : bins         ,
            "method"       : "LBL"        ,
            }

        return None

def _spec_shift(
    X      : tuple[np.ndarray],
    *param : tuple[float]
    ) -> np.ndarray:
    """Spectrum shifted.

    Parameters
    ----------
    X : tuple[np.ndarray]
        wavelength, flux and gradient arrays
    *param : tuple[float]
        free parameter(s): z (RV divided by speed of light)

    Returns
    -------
    np.ndarray
        shifted spectrum
    """

    wave, flux, grad = X
    z = param
    S = flux - grad*wave*z

    return S

def _spec_shift_scale(
    X      : tuple[np.ndarray],
    *param : tuple[float]
    ) -> np.ndarray:
    """Spectrum shifted and scaled.

    Parameters
    ----------
    X : tuple[np.ndarray]
        wavelength, flux and gradient arrays
    *param : tuple[float]
        free parameter(s): z (RV divided by speed of light) and A (scaling factor)

    Returns
    -------
    np.ndarray
        shifted and scaled spectrum
    """

    wave, flux, grad = X
    z, A = param
    S = (flux - grad*wave*z)*A

    return S