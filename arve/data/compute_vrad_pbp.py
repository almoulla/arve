from   astropy.stats     import sigma_clip
import numpy             as     np
import pandas            as     pd
from   scipy.interpolate import interp1d
from   scipy.optimize    import curve_fit
from   tqdm              import tqdm
import warnings
warnings.filterwarnings("ignore")

class compute_vrad_pbp:

    def compute_vrad_pbp(self, criteria:list=None, exclude_tellurics=True, scale=True, Niter=1, bins:"list[list]"=None) -> None:
        """Compute radial velocities (RVs) from spectral data.

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
        Nspec                        =  self.spec["Nspec"]
        Nord                         =  self.spec["Nord"]

        # master spectrum gradient
        mast_grad_val = np.zeros_like(mast_wave_val)
        for i in range(Nord):
            mast_grad_val[i] = np.gradient(mast_flux_val[i])/np.gradient(mast_wave_val[i])
        
        # flux and gradient functions
        flux_func = np.zeros(Nord, dtype=object)
        grad_func = np.zeros(Nord, dtype=object)
        for i in range(Nord):
            flux_func[i] = interp1d(mast_wave_val[i], mast_flux_val[i], kind="cubic")
            grad_func[i] = interp1d(mast_wave_val[i], mast_grad_val[i], kind="cubic")

        # read constants
        c = self.arve.functions.constants["c"]

        # read mask
        mask = self.aux_data["mask"]

        # central, left and right wavelengths
        wc = np.zeros(Nord, dtype=object)
        wl = np.zeros(Nord, dtype=object)
        wr = np.zeros(Nord, dtype=object)
        for i in range(Nord):
            wc[i] = np.array(mask[i]["wave"  ])
            wl[i] = np.array(mask[i]["wave_l"])
            wr[i] = np.array(mask[i]["wave_r"])

        # temperature bins
        if bins is None:
            bins = [[np.min(temp),np.max(temp)]]

        # nr. of lines and temperature bins
        Nline = np.zeros(Nord, dtype=object)
        for i in range(Nord):
            Nline[i] = len(wc[i])
        Ntemp = len(bins)

        # exclude tellurics
        if exclude_tellurics:
            if criteria is None:
                criteria = ["tell"]
            else:
                criteria.append("tell")

        # lines which satisfy criteria
        idx_crit = np.zeros(Nord, dtype=object)
        for i in range(Nord):
            idx_crit[i] = np.ones_like(wc[i], dtype=bool)
            if criteria is not None:
                for j in range(len(criteria)):
                    crit = np.array(mask[i]["crit_"+criteria[j]])
                    idx_crit[i] *= crit

        # index of lines
        idx_line = np.zeros(Nord, dtype=object)
        for i in range(Nord):
            idx_line[i] = np.zeros(Nline[i], dtype=object)
            for j in range(Nline[i]):
                il             = np.searchsorted(wave_val[i], wl[i][j], side="right")
                ir             = np.searchsorted(wave_val[i], wr[i][j], side="right")
                idx_line[i][j] = np.arange(il, ir)

        # index of temperature bins
        idx_temp = np.zeros(Nord, dtype=object)
        for i in range(Nord):
            idx_temp[i] = np.zeros(Ntemp, dtype=object)
            for j in range(Ntemp):
                il             = temp[i]>=bins[j][0]
                ir             = temp[i]<=bins[j][1]
                idx_temp[i][j] = np.where(il & ir)[0]

        # index of temperature-segmented line parts
        idx = np.zeros(Nord, dtype=object)
        for i in range(Nord):
            idx[i] = np.zeros((Nline[i],Ntemp), dtype=object)
            for j in range(Nline[i]):
                for k in range(Ntemp):
                    idx[i][j,k] = np.intersect1d(idx_line[i][j], idx_temp[i][k])

        # empty arrays for RV values and errors
        vrad_val_arr = np.zeros(Nord, dtype=object)
        vrad_err_arr = np.zeros(Nord, dtype=object)
        for i in range(Nord):
            vrad_val_arr[i] = np.zeros((Nline[i],Ntemp,Nspec))*np.nan
            vrad_err_arr[i] = np.zeros((Nline[i],Ntemp,Nspec))*np.nan

        # loop spectra
        print("Analyzed spectra:")
        for i in tqdm(range(Nspec)):

            # read spectrum
            _, flux_val, flux_err = self.read_spec(i)

            # loop orders
            for j in range(Nord):

                # loop lines
                for k in range(Nline[j]):

                    # check if line satisfies criteria
                    if idx_crit[j][k] == True:

                        # loop temperature bins
                        for l in range(Ntemp):

                            # initial guess
                            vrad_val = 0

                            # wavelength, flux and flux uncertainty for line
                            wave_spec_line = wave_val[j,idx[j][k,l]]
                            flux_spec_line = flux_val[j,idx[j][k,l]]
                            ferr_spec_line = flux_err[j,idx[j][k,l]]

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
                                    flux_mast_inter = flux_func[j](wave_mast_inter)
                                    grad_mast_inter = grad_func[j](wave_mast_inter)

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
                                vrad_val_arr[j][k,l,i] = vrad_val
                                vrad_err_arr[j][k,l,i] = vrad_err
                            
                            # if unsuccessful, continue
                            except:

                                continue

        # outlier rejection
        for i in range(Nspec):
            *_, clip_val_min, clip_val_max = sigma_clip(np.concatenate(vrad_val_arr)[:,:,i], maxiters=None, return_bounds=True)
            *_, clip_err_min, clip_err_max = sigma_clip(np.concatenate(vrad_err_arr)[:,:,i], maxiters=None, return_bounds=True)
            for j in range(Nord):
                for k in range(Nline[j]):
                    for l in range(Ntemp):
                        if (vrad_val_arr[j][k,l,i] < clip_val_min) | (vrad_val_arr[j][k,l,i] > clip_val_max):
                            vrad_val_arr[j][k,l,i] = np.nan
                        if (vrad_err_arr[j][k,l,i] < clip_err_min) | (vrad_err_arr[j][k,l,i] > clip_err_max):
                            vrad_err_arr[j][k,l,i] = np.nan

        # weighted average of all valid lines per order
        vrad_val_ord = np.zeros((Nspec,Nord,Ntemp))*np.nan
        vrad_err_ord = np.zeros((Nspec,Nord,Ntemp))*np.nan
        for i in range(Nspec):
            for j in range(Nord):
                for k in range(Ntemp):
                    idx = ~np.isnan(vrad_val_arr[j][:,k,i]) & (np.abs(vrad_err_arr[j][:,k,i])>0)
                    if np.sum(idx) > 0:
                        vrad_val_ord[i,j,k] = np.average(vrad_val_arr[j][idx,k,i], weights=1/vrad_err_arr[j][idx,k,i]**2)
                        vrad_err_ord[i,j,k] = np.sqrt(1/np.sum(1/vrad_err_arr[j][idx,k,i]**2))

        # weighted average of all orders
        vrad_val = np.zeros((Nspec,Ntemp))*np.nan
        vrad_err = np.zeros((Nspec,Ntemp))*np.nan
        for i in range(Nspec):
            for j in range(Ntemp):
                idx = ~np.isnan(vrad_val_ord[i,:,j]) & (np.abs(vrad_err_ord[i,:,j])>0)
                if np.sum(idx) > 0:
                    vrad_val[i,j] = np.average(vrad_val_ord[i,idx,j], weights=1/vrad_err_ord[i,idx,j]**2)
                    vrad_err[i,j] = np.sqrt(1/np.sum(1/vrad_err_ord[i,idx,j]**2))

        # save RV data
        self.vrad = {
            "vrad_val"    : vrad_val    ,
            "vrad_err"    : vrad_err    ,
            "vrad_val_ord": vrad_val_ord,
            "vrad_err_ord": vrad_err_ord,
            "vrad_val_pbp": vrad_val_arr,
            "vrad_err_pbp": vrad_err_arr,
            "bins"        : bins        ,
            "method"      : "PBP"       ,
            }

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