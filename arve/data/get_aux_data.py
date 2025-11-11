import importlib.resources
import os
import urllib.request

import numpy                 as     np
import pandas                as     pd
from   scipy.interpolate     import interp1d
from   scipy.ndimage.filters import convolve

from typing import Literal

class get_aux_data:

    def get_aux_data(
        self,
        mask_path       : str                   | None = None ,
        mask_medium     : Literal["vac", "air"]        = "vac",
        tell_wave       : np.ndarray            | None = None ,
        tell_flux       : np.ndarray            | None = None ,
        tell_lim        : float                        = 0.99 ,
        exclude_regions : list[list[float]]     | None = None ,
        ) -> None:
        """Get auxiliary data.

        Parameters
        ----------
        mask_path : str | None, optional
            path to self-provided line mask (must be a CSV file where the relevant columns are named the same as in the package-provided masks), by default None
        mask_medium : Literal["vac", "air"], optional
            medium of mask wavelengths, by default "vac"
        tell_wave : np.ndarray | None, optional
            wavelength values of normalized telluric spectrum (to replace the general model), by default None
        tell_flux : np.ndarray | None, optional
            flux values of normalized telluric spectrum (to replace the general model), by default None
        tell_lim : float, optional
            normalized flux limit on telluric features (to add a criterion column for stellar lines unaffected by tellurics), by default 0.99
        exclude_regions : list[list[float]] | None, optional
            wavelength intervals to be excluded (to add a criterion column for stellar lines outside the excluded regions), by default None

        Returns
        -------
        None
            None
        """
        
        # read wavelengh values
        wave_val = self.spec["wave_val"]

        # read resolution and medium
        resolution = self.spec["resolution"]
        medium     = self.spec["medium"]

        # read systematic velocity and maximum BERV
        vrad_sys = self.arve.star.stellar_parameters["vrad_sys"]
        berv_max = self.arve.star.stellar_parameters["berv_max"]

        # paths to auxiliary data
        resource_aux_data = importlib.resources.files("arve.aux_data")
        with importlib.resources.as_file(resource_aux_data) as path:
            path_aux_data = str(path) + "/"
        path_aux_mask = path_aux_data+"masks/"
        path_aux_spec = path_aux_data+"spectra/"
        path_aux_tell = path_aux_data+"tellurics/"
        path_aux_wave = path_aux_data+"wavelengths/"

        # search file names and spectral types among masks
        files             = sorted(os.listdir(path_aux_mask))
        files             = [file for file in files if file.endswith(".csv.zip")]
        sptypes_data      = [file.split(".")[0] for file in files]
        sptypes_data_num  = np.array([self.arve.functions.sptype_to_num(sptype=sptype) for sptype in sptypes_data])

        # spectral type as number
        sptype_star     = self.arve.star.stellar_parameters["sptype"]
        sptype_star_num = self.arve.functions.sptype_to_num(sptype=sptype_star)

        # file of closest spectral type available in auxiliary data
        idx_file = np.argmin(np.abs(sptype_star_num-sptypes_data_num))
        file     = files[idx_file]

        # read closest mask
        if mask_path is None:
            mask      = pd.read_csv(path_aux_mask+file)
            mask_name = file
        # read self-provided mask
        else:
            mask      = pd.read_csv(mask_path)
            mask_name = mask_path.split("/")[-1]
        
        # convert wavelengths to correct medium
        if (medium=="air") & (mask_medium=="vac"):
            mask["wave"  ] = self.arve.functions.convert_vac_to_air(mask["wave"  ])
            mask["wave_l"] = self.arve.functions.convert_vac_to_air(mask["wave_l"])
            mask["wave_u"] = self.arve.functions.convert_vac_to_air(mask["wave_u"])
        if (medium=="vac") & (mask_medium=="air"):
            mask["wave"  ] = self.arve.functions.convert_air_to_vac(mask["wave"  ])
            mask["wave_l"] = self.arve.functions.convert_air_to_vac(mask["wave_l"])
            mask["wave_u"] = self.arve.functions.convert_air_to_vac(mask["wave_u"])

        # keep only mask lines within wavelength range
        mask_dict = {}
        for i in range(self.spec["N_ord"]):
            idx_wave_l   = np.searchsorted(mask["wave_l"], wave_val[i][ 0])
            idx_wave_u   = np.searchsorted(mask["wave_u"], wave_val[i][-1])
            idx_wave     = np.arange(idx_wave_l, idx_wave_u)
            mask_dict[i] = mask.iloc[idx_wave].reset_index(drop=True)

        # download closest spectrum from GitHub if it does not already exist in the package directory
        if not os.path.exists(path_aux_spec+file):
            path_github = "https://raw.githubusercontent.com/almoulla/arve/main/arve/aux_data/spectra/"
            urllib.request.urlretrieve(path_github+file, path_aux_spec+file)

        # read closest spectrum
        wave = pd.read_csv(path_aux_wave+"WAVE.csv.zip")["wave"]
        spec = pd.read_csv(path_aux_spec+file)
        if medium == "air":
            wave = self.arve.functions.convert_vac_to_air(wave)
        spec.insert(0, "wave", wave)

        # convolve spectrum
        if resolution is not None:
            spec["flux"] = _convolve_gaussian(spec["wave"].values, spec["flux"].values, resolution)

        # interpolate spectrum
        spec_wave = wave_val
        spec_flux = np.array([interp1d(spec["wave"], spec["flux"], kind="cubic", bounds_error=False)(spec_wave[i]) for i in range(self.spec["N_ord"])])
        spec_temp = np.array([interp1d(spec["wave"], spec["temp"], kind="cubic", bounds_error=False)(spec_wave[i]) for i in range(self.spec["N_ord"])])
        spec_dict = {
        "wave": spec_wave,
        "flux": spec_flux,
        "temp": spec_temp
        }

        # find lower and upper index of lines
        for i in range(self.spec["N_ord"]):
            if "wave" in mask_dict[0].columns:
                col_ref = "wave"
            if "flux" in mask_dict[0].columns:
                col_ref = "flux"
            mask_dict[i].insert(mask_dict[i].columns.get_loc(col_ref+"_l")+2, "idx_l", np.searchsorted(spec_dict["wave"][i], mask_dict[i]["wave_l"])  )
            mask_dict[i].insert(mask_dict[i].columns.get_loc(col_ref+"_u")+2, "idx_u", np.searchsorted(spec_dict["wave"][i], mask_dict[i]["wave_u"])-1)
            mask_dict[i]["wave_l"] = spec_dict["wave"][i][mask_dict[i]["idx_l"]]
            mask_dict[i]["wave_u"] = spec_dict["wave"][i][mask_dict[i]["idx_u"]]
        
        # construct telluric spectrum with a general model interpolated and convolved as the stellar spectrum
        if (tell_wave is None) & (tell_flux is None):

            # read telluric spectrum
            wave = pd.read_csv(path_aux_wave+"WAVE.csv.zip")["wave"]
            tell = pd.read_csv(path_aux_tell+"TELL.csv.zip")
            if medium == "air":
                wave = self.arve.functions.convert_vac_to_air(wave)
            wave = self.arve.functions.doppler_shift(wave=wave, v=-vrad_sys)
            tell.insert(0, "wave", wave)

            # convolve telluric spectrum
            if resolution is not None:
                tell["flux"] = _convolve_gaussian(tell["wave"].values, tell["flux"].values, resolution)

            # interpolate telluric spectrum
            tell_wave = wave_val
            tell_flux = np.array([interp1d(tell["wave"], tell["flux"], kind="cubic", bounds_error=False)(tell_wave[i]) for i in range(self.spec["N_ord"])])
            tell_dict = {
            "wave": tell_wave,
            "flux": tell_flux
            }

        # or use provided telluric spectrum (assumed to have the same resolution as the stellar spectrum but not necessarily sampled on the same wavelength grid)
        else:

            # concatenate stellar spectrum
            tell_wave_1d = np.concatenate(tell_wave)
            tell_flux_1d = np.concatenate(tell_flux)
            idx_sort = np.argsort(tell_wave_1d)

            # create DataFrame with flattened telluric spectrum
            tell = pd.DataFrame()
            tell["wave"] = tell_wave_1d[idx_sort]
            tell["flux"] = tell_flux_1d[idx_sort]

            # create dictionary with interpolated telluric spectrum
            tell_wave = wave_val
            tell_flux = np.array([interp1d(tell_wave[i], tell_flux[i], kind="cubic", bounds_error=False)(tell_wave[i]) for i in range(self.spec["N_ord"])])
            tell_dict = {
            "wave": tell_wave,
            "flux": tell_flux
            }

        # find lower and upper bounds of telluric bands
        tell_wave = np.array(tell["wave"])
        tell_flux = np.array(tell["flux"])
        idx_above = np.where(tell_flux>=tell_lim)[0]
        idx_below = np.where(tell_flux< tell_lim)[0]
        if len(idx_above) == 0:
            tell_wave_l = np.array([tell_wave[ 0]])
            tell_wave_u = np.array([tell_wave[-1]])
        elif len(idx_below) == 0:
            tell_wave_l = np.array([])
            tell_wave_u = np.array([])
        else:
            idx_l = idx_above[:-1][np.diff(idx_above)>1]
            idx_u = idx_below[:-1][np.diff(idx_below)>1]
            if tell_flux[ 0] <  tell_lim:
                idx_l = np.append(0, idx_l)
            if tell_flux[-1] >= tell_lim:
                idx_u = np.append(idx_u, idx_below[-1])
            if tell_flux[-1] <  tell_lim:
                idx_l = np.append(idx_l, idx_above[-1])
                idx_u = np.append(idx_u, len(tell_flux)-1)
            tell_wave_l = self.arve.functions.doppler_shift(wave=tell_wave[idx_l], v=-berv_max)
            tell_wave_u = self.arve.functions.doppler_shift(wave=tell_wave[idx_u], v= berv_max)

        # create new telluric DataFrame with lower and upper bounds of telluric bands
        band = pd.DataFrame()
        band["wave_l"] = tell_wave_l
        band["wave_u"] = tell_wave_u

        # keep telluric bands which overlap with data
        band_dict = {}
        for i in range(self.spec["N_ord"]):
            idx_wave     = (band["wave_l"]>wave_val[i][-1]) | (band["wave_u"]<wave_val[i][0]) == False
            band_dict[i] = band[idx_wave].reset_index(drop=True)
            band_dict[i]["wave_l"][band_dict[i]["wave_l"]<wave_val[i][ 0]] = wave_val[i][ 0]
            band_dict[i]["wave_u"][band_dict[i]["wave_u"]>wave_val[i][-1]] = wave_val[i][-1]

        # create new telluric DataFrame with lower and upper bounds of excluded regions
        if exclude_regions is not None:
            exclude_regions = np.array(exclude_regions)
            excl = pd.DataFrame()
            excl["wave_l"] = exclude_regions[:,0]
            excl["wave_u"] = exclude_regions[:,1]

        # keep excluded regions which overlap with data
        excl_dict = {}
        for i in range(self.spec["N_ord"]):
            if exclude_regions is not None:
                idx_wave     = (excl["wave_l"]>wave_val[i][-1]) | (excl["wave_u"]<wave_val[i][0]) == False
                excl_dict[i] = excl[idx_wave].reset_index(drop=True)
                excl_dict[i]["wave_l"][excl_dict[i]["wave_l"]<wave_val[i][ 0]] = wave_val[i][ 0]
                excl_dict[i]["wave_u"][excl_dict[i]["wave_u"]>wave_val[i][-1]] = wave_val[i][-1]
            else:
                excl_dict[i] = pd.DataFrame(columns=["wave_l", "wave_u"])

        # criterion: lines not within telluric bands
        for i in range(self.spec["N_ord"]):
            band_wave_l = band_dict[i]["wave_l"].values
            band_wave_u = band_dict[i]["wave_u"].values
            mask_wave_l = mask_dict[i]["wave_l"].values
            mask_wave_u = mask_dict[i]["wave_u"].values
            mask_dict[i]["crit_tell"] = np.sum([(mask_wave_l>band_wave_u[j]) | (mask_wave_u<band_wave_l[j]) for j in range(len(band_dict[i]))], axis=0) == len(band_dict[i])

        # criterion: lines not within excluded regions
        for i in range(self.spec["N_ord"]):
            excl_wave_l = excl_dict[i]["wave_l"].values
            excl_wave_u = excl_dict[i]["wave_u"].values
            mask_wave_l = mask_dict[i]["wave_l"].values
            mask_wave_u = mask_dict[i]["wave_u"].values
            mask_dict[i]["crit_excl"] = np.sum([(mask_wave_l>excl_wave_u[j]) | (mask_wave_u<excl_wave_l[j]) for j in range(len(excl_dict[i]))], axis=0) == len(excl_dict[i])

        # save
        self.aux_data = {"name": mask_name,
                         "mask": mask_dict,
                         "spec": spec_dict,
                         "tell": tell_dict,
                         "band": band_dict,
                         "excl": excl_dict
                        }

        return None

def _convolve_gaussian(
    wave       : np.ndarray,
    flux       : np.ndarray,
    resolution : float
    ) -> np.ndarray:
    """Convolve spectrum with a Gaussian instrumental profile. Adapted from: https://pysme-astro.readthedocs.io/en/latest/_modules/pysme/broadening.html

    Parameters
    ----------
    wave : np.ndarray
        wavelength values
    flux : np.ndarray
        flux values
    resolution : float
        instrumental resolution

    Returns
    -------
    np.ndarray
        flux values convolved with a Gaussian instrumental profile
    """

    # nr. of wavelength points
    Nwave = len(wave)

    # half-width at half-maximum
    hwhm = 0.5*wave[0]/resolution

    # uniform dispersion
    wave_uan = wave[-1]-wave[0]
    wave_del = wave_uan/(Nwave-1)
    
    # nr. of points in half Gaussian
    Nhalf = int(4/np.sqrt(2*np.log(2))*hwhm/wave_del)
    
    # nr. of points in fill Gaussian (odd)
    Nfull = 2*Nhalf+1
    
    # wavelength grid of Gaussian
    wave_grid = wave_del*(np.arange(Nfull)-(Nfull-1)/2)
    
    # standard deviation
    sigma = hwhm/np.sqrt(2*np.log(2))
    
    # normalized Gaussian
    gauss  = np.exp(-(wave_grid/sigma)**2/2)
    gauss /= np.sum(gauss)

    # convolve spectrum
    flux = convolve(flux, gauss, mode="nearest")

    return flux