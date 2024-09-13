import numpy             as     np
import pandas            as     pd
from   scipy.interpolate import interp1d
from   scipy.optimize    import curve_fit
from   tqdm              import tqdm

class compute_vrad_ccf:

    def compute_vrad_ccf(self, weight:str=None, criteria:list=None, exclude_tellurics=True, exclude_regions:bool=True, vgrid:list=None) -> None:
        """Compute radial velocities (RVs) from spectral data.

        :param weight: column name of weight, defaults to None
        :type weight: str, optional
        :param criteria: criteria to apply (must be columns with prefix "crit_"), defaults to None
        :type criteria: list, optional
        :param exclude_tellurics: exclude telluric bands, defaults to True
        :type exclude_tellurics: bool, optional
        :param exclude_regions: exclude wavelength intervals, defaults to True
        :type exclude_regions: bool, optional
        :param vgrid: velocity grid, in the format [start,stop,step] and in units of km/s, on which to evaluate the CCF, defaults to None
        :type vgrid: list, optional
        :return: None
        :rtype: None
        """

        # read data
        wave_val = self.spec["wave_val"]
        Nord     = self.spec["Nord"]

        # read constants
        c = self.arve.functions.constants["c"]

        # read mask
        mask = self.aux_data["mask"]

        # central wavelength
        wc = np.zeros(Nord, dtype=object)
        for i in range(Nord):
            wc[i] = np.array(mask[i]["wave"])

        # weights
        w = np.zeros(Nord, dtype=object)
        for i in range(Nord):
            if weight is None:
                w[i] = np.ones(len(mask[i]))
            else:
                w[i] = np.array(mask[i][weight])

        # exclude tellurics
        if exclude_tellurics:
            if criteria is None:
                criteria = ["tell"]
            else:
                criteria.append("tell")
        
        # exclude regions
        if exclude_regions:
            if criteria is None:
                criteria = ["excl"]
            else:
                criteria.append("excl")

        # lines which satisfy criteria
        if criteria is not None:
            for i in range(Nord):
                idx = np.ones_like(wc[i], dtype=bool)
                for j in range(len(criteria)):
                    crit = np.array(mask[i]["crit_"+criteria[j]])
                    idx *= crit
                wc[i] = wc[i][idx]
                w [i] = w [i][idx]

        # RV shifts
        if vgrid is None:
            wave_val_flat = np.concatenate(wave_val)
            vrads_start   = -20
            vrads_stop    =  20
            vrads_step    = np.nanmedian(np.diff(wave_val_flat)/wave_val_flat[:-1])*c
            vrads         = np.arange(vrads_start, vrads_stop, vrads_step)
        else:
            vrads         = np.arange(vgrid[0],vgrid[1]+vgrid[2]/2,vgrid[2])

        # keep mask lines within spectrum overlap
        for i in range(Nord):
            idx = (self.arve.functions.doppler_shift(wave=wc[i], v=np.min(vrads)) > np.min(wave_val[i])) & \
                  (self.arve.functions.doppler_shift(wave=wc[i], v=np.max(vrads)) < np.max(wave_val[i]))
            wc[i] = wc[i][idx]
            w [i] = w [i][idx]

        # normalize weights
        for i in range(Nord):
            w[i] = w[i]/np.sum(w[i])
        
        # nr. of spectra, RV shifts and lines
        if self.spec["path"] is None:
            Nspec = len(self.spec["flux_val"])
        else:
            Nspec = len(self.spec["files"])
        Nvrad = len(vrads)

        # empty arrays for indices
        ir = np.zeros((Nord,Nvrad), dtype=object)
        il = np.zeros((Nord,Nvrad), dtype=object)

        # empty arrays for fractions
        fr = np.zeros((Nord,Nvrad), dtype=object)
        fl = np.zeros((Nord,Nvrad), dtype=object)

        # loop spectral orders
        for i in range(Nord):

            # loop RV shifts
            for j in range(Nvrad):

                # shifted central wavelengths
                wc_shift = self.arve.functions.doppler_shift(wave=wc[i], v=vrads[j])

                # right and left indices
                ir[i,j] = np.searchsorted(wave_val[i], wc_shift, side="right")
                il[i,j] = ir[i,j] - 1

                # right and left fractions
                fr[i,j] = (wc_shift-wave_val[i][il[i,j]])/(wave_val[i][ir[i,j]]-wave_val[i][il[i,j]])
                fl[i,j] = 1 - fr[i,j]

        # empty arrays for CCF values and errors
        ccf_val = np.zeros((Nspec,Nord+1,Nvrad))
        ccf_err = np.zeros((Nspec,Nord+1,Nvrad))

        # empty arrays for CCF moments
        vrad_val_arr = np.zeros((Nspec,Nord+1))*np.nan
        vrad_err_arr = np.zeros((Nspec,Nord+1))*np.nan
        cont_val_arr = np.zeros((Nspec,Nord+1))*np.nan
        cont_err_arr = np.zeros((Nspec,Nord+1))*np.nan
        ampl_val_arr = np.zeros((Nspec,Nord+1))*np.nan
        ampl_err_arr = np.zeros((Nspec,Nord+1))*np.nan
        fwhm_val_arr = np.zeros((Nspec,Nord+1))*np.nan
        fwhm_err_arr = np.zeros((Nspec,Nord+1))*np.nan

        # loop spectra
        print("Analyzed spectra:")
        for i in tqdm(range(Nspec)):

            # read spectrum
            _, flux_val, flux_err = self.read_spec(i)

            # loop orders
            for j in range(Nord+1):

                try:

                    # single order
                    if j < Nord:

                        # loop RV shifts
                        for k in range(Nvrad):

                            # CCF value and error
                            ccf_val[i,j,k] = np.nansum(( flux_val[j][il[j,k]]*fl[j,k]    + flux_val[j][ir[j,k]]*fr[j,k]    )*w[j]   )
                            ccf_err[i,j,k] = np.nansum(((flux_err[j][il[j,k]]*fl[j,k])**2+(flux_err[j][ir[j,k]]*fr[j,k])**2)*w[j]**2)**(1/2)

                    # sum of all orders
                    else:

                        # order-summed CCF values and errors
                        ccf_val[i,j] = np.nansum(ccf_val[i]   , axis=0)
                        ccf_err[i,j] = np.nansum(ccf_err[i]**2, axis=0)**(1/2)

                    # initial guess on Gaussian parameters
                    i_min = np.argmin(ccf_val[i,j])
                    i_max = np.argmax(ccf_val[i,j])
                    C0 = ccf_val[i,j,i_max]
                    a0 = 1-ccf_val[i,j,i_min]/C0
                    b0 = vrads[i_min]
                    c0 = (b0-interp1d(ccf_val[i,j,:i_min], vrads[:i_min], kind="cubic")((ccf_val[i,j,i_min]+ccf_val[i,j,i_max])/2))*2
                    p0 = (C0, a0, b0, c0)

                    # CCF parameters with fitted Gaussian
                    param_val, param_cov = curve_fit(self.arve.functions.inverted_gaussian, vrads, ccf_val[i,j], sigma=ccf_err[i,j], p0=p0)
                    param_err            = np.sqrt(np.diag(param_cov))

                    # RV value and error
                    vrad_val_arr[i,j] = param_val[2]
                    vrad_err_arr[i,j] = 1/np.sqrt(np.sum(1/(ccf_err[i,j]*np.gradient(vrads)/np.gradient(ccf_val[i,j]))**2))

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
        vrad_val = np.zeros(Nspec)*np.nan
        vrad_err = np.zeros(Nspec)*np.nan
        cont_val = np.zeros(Nspec)*np.nan
        cont_err = np.zeros(Nspec)*np.nan
        ampl_val = np.zeros(Nspec)*np.nan
        ampl_err = np.zeros(Nspec)*np.nan
        fwhm_val = np.zeros(Nspec)*np.nan
        fwhm_err = np.zeros(Nspec)*np.nan

        # loop spectra
        for i in range(Nspec):

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
        bisector_flux_val = np.arange(0.00, 1.00, 0.01)
        Nbis              = len(bisector_flux_val)

        # emtpty arrays for bisector and BIS value and error
        bisector_vrad_val = np.zeros((Nspec,Nbis))*np.nan
        bisector_vrad_err = np.zeros((Nspec,Nbis))*np.nan
        bis_val           = np.zeros( Nspec      )*np.nan
        bis_err           = np.zeros( Nspec      )*np.nan

        # loop spectra
        for i in range(Nspec):

            try:

                # CCF RV error
                ccf_err_vrad = (ccf_err[i,-1]*np.gradient(vrads)/np.gradient(ccf_val[i,-1]))**2

                # bisector value
                bisector_vrad_val_l  = interp1d(ccf_val[i,-1][vrads<vrad_val[i]]/cont_val[i], vrads[vrads<vrad_val[i]], bounds_error=False)(bisector_flux_val)
                bisector_vrad_val_r  = interp1d(ccf_val[i,-1][vrads>vrad_val[i]]/cont_val[i], vrads[vrads>vrad_val[i]], bounds_error=False)(bisector_flux_val)
                bisector_vrad_val[i] = (bisector_vrad_val_l+bisector_vrad_val_r)/2

                # bisector error
                bisector_vrad_err_l  = interp1d(vrads[vrads<vrad_val[i]], ccf_err_vrad[vrads<vrad_val[i]], bounds_error=False)(bisector_vrad_val_l)
                bisector_vrad_err_r  = interp1d(vrads[vrads>vrad_val[i]], ccf_err_vrad[vrads>vrad_val[i]], bounds_error=False)(bisector_vrad_val_r)
                bisector_vrad_err[i] = np.sqrt(bisector_vrad_err_l**2 + bisector_vrad_err_r**2)/2

                # bisector normalized flux
                bisector_flux_norm = 1-(1-bisector_flux_val)/ampl_val[i]

                # bisector index for bottom and top
                idx_bisector_bot = (bisector_flux_norm>0.1) & (bisector_flux_norm<0.4)
                idx_bisector_top = (bisector_flux_norm>0.6) & (bisector_flux_norm<0.9)

                # BIS bottom and top value and error
                bis_bot_val = np.average(bisector_vrad_val[i][idx_bisector_bot])
                bis_bot_err = np.average(bisector_vrad_err[i][idx_bisector_bot])/np.sqrt(np.sum(idx_bisector_bot))
                bis_top_val = np.average(bisector_vrad_val[i][idx_bisector_top])
                bis_top_err = np.average(bisector_vrad_err[i][idx_bisector_top])/np.sqrt(np.sum(idx_bisector_top))

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
            "vrads"            : vrads            ,
            "ccf_val"          : ccf_val          ,
            "ccf_err"          : ccf_err          ,
            "bisector_vrad_val": bisector_vrad_val,
            "bisector_vrad_err": bisector_vrad_err,
            "bisector_flux_val": bisector_flux_val,
            "cont_val"         : cont_val         ,
            "cont_err"         : cont_err         ,
            "ampl_val"         : ampl_val         ,
            "ampl_err"         : ampl_err         ,
            "fwhm_val"         : fwhm_val         ,
            "fwhm_err"         : fwhm_err         ,
            "bis_val"          : bis_val          ,
            "bis_err"          : bis_err          ,
            }

        return None

def _weighted_average(val_arr, err_arr):

    idx = ~np.isnan(val_arr) & ~np.isnan(err_arr) & (err_arr>0)
    if np.sum(idx) > 0:
        val_avg = np.average(val_arr[idx], weights=1/err_arr[idx]**2)
        err_avg = np.sqrt(1/np.sum(1/err_arr[idx]**2))
    else:
        val_avg = np.nan
        err_avg = np.nan

    return val_avg, err_avg