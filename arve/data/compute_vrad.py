from   numba          import njit
import numpy          as     np
import os
import pkg_resources
from   scipy.optimize import curve_fit

def compute_vrad(self, methods:list=["CCF"]) -> None:
    """Compute radial velocities (RVs) from spectral data.

    :param methods: extraction methods, defaults to ["CCF"]
    :type methods: list, optional
    :return: None
    :rtype: None
    """

    # read spectral data and units
    time, wave, flux_val, flux_err = [self.spec[key] for key in ["time", "wave", "flux_val", "flux_err"]]
    time_unit, = [self.spec[key] for key in ["time_unit"]]

    # nr. of spectra
    Nspec = flux_val.shape[0]

    # method CCF
    if "CCF" in methods:
    
        # RV shifts
        vrads = np.arange(-40,40,0.25)
        Nvrad = len(vrads)

        # search masks
        path_aux_data = pkg_resources.resource_filename("arve", "aux_data/")
        masks         = os.listdir(path_aux_data+"masks/")
        masks         = [mask for mask in masks if mask.endswith(".mask")]
        sptype_masks  = [mask.split('_')[0] for mask in masks]

        # convert spectral types to numbers
        sptype_num       =           self.arve.functions.sptype_to_num(sptype=self.arve.star.stellar_parameters["sptype"])
        sptype_num_masks = np.array([self.arve.functions.sptype_to_num(sptype=sptype) for sptype in sptype_masks])

        # read closest mask
        idx_mask = np.argmin(np.abs(sptype_num-sptype_num_masks))
        mask     = np.genfromtxt(path_aux_data+"masks/"+masks[idx_mask])

        # central wavelengths (wc) and weights (w)
        wc = mask[:,0]
        w  = mask[:,1]

        # wavelength overlap between mask and shifted spectrum
        wc_min = max(min(wc), self.arve.functions.doppler_shift(wave=min(wave), v=min(vrads)))
        wc_max = min(max(wc), self.arve.functions.doppler_shift(wave=max(wave), v=max(vrads)))
        
        # keep mask lines within overlap and re-normalize weights
        idx = (wc > wc_min) & (wc < wc_max)
        wc  = wc[idx]
        w   = w [idx]
        w   = w/np.sum(w)

        # empty arrays for RV values and errors
        vrad_val = np.zeros(Nspec)
        vrad_err = np.zeros(Nspec)

        # loop spectra
        for i in range(Nspec):

            # empty arrays for CCF values and errors
            ccf_val = np.zeros(Nvrad)
            ccf_err = np.zeros(Nvrad)

            # loop RV shifts
            for j in range(Nvrad):

                # RV shift
                v = vrads[j]

                # shifted wavelengths and interpolated spectra
                wave_shift     = self.arve.functions.doppler_shift(wave=wave, v=v)
                flux_inter_val = np.interp(wc, wave_shift, flux_val[i])
                flux_inter_err = np.interp(wc, wave_shift, flux_err[i])

                # CCF values and errors
                ccf_val[j], ccf_err[j] = _func_ccf(flux_inter_val, flux_inter_err, w)

            # initial guess
            C0 = np.max(ccf_val)
            a0 = C0-np.min(ccf_val)
            b0 = vrads[np.argmin(ccf_val)]
            c0 = (b0-vrads[np.where(ccf_val<(C0-a0/2))[0][0]])*2
            p0 = (C0, a0, b0, c0)

            # CCF moments with fitted Gaussian
            ccf_mom, _ = curve_fit(_func_inverted_gauss, vrads, ccf_val, sigma=ccf_err, p0=p0)

            # RV value and error
            vrad_val[i] = ccf_mom[2]
            vrad_err[i] = 1/np.sqrt(np.sum(1/np.abs(ccf_err*np.gradient(vrads)/np.gradient(ccf_val))**2))

        # save RV data
        self.vrad = {"time": time, "vrad_val": vrad_val, "vrad_err": vrad_err, "time_unit": time_unit, "vrad_unit": "km/s", "mask": masks[idx_mask][:2]}

    return None

@njit()
def _func_ccf(flux_val:list, flux_err:list, w:list) -> tuple:
    """Cross-correlation function (CCF).

    :param flux_val: flux values
    :type flux_val: list
    :param flux_err: flux errors
    :type flux_err: list
    :param w: weights
    :type w: list
    :return: CCF values and errors
    :rtype: tuple
    """

    # CCF values and errors
    ccf_val = np.sum(flux_val*w)
    ccf_err = np.sum(flux_err**2*w)**(1/2)

    return ccf_val, ccf_err

def _func_inverted_gauss(x:list, *params:tuple) -> list:
    """Inverted Gaussian.

    :param x: RV array
    :type x: list
    :param params: tuple with continuum, contrast, RV and FWHM
    :type params: tuple of floats
    :return: inverted Gaussian evaluated at x
    :rtype: list
    """
    
    # unpack parameters
    continuum, contrast, RV, FWHM = params
    
    # rename parameters
    C = continuum
    a = contrast
    b = RV
    c = FWHM
    
    # scale FWHM into sigma
    c /= 2*np.sqrt(2*np.log(2))

    return C - a*np.exp(-((x-b)/c)**2/2)