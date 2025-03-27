import importlib.resources
import os
import urllib.request

import numpy                 as     np
import pandas                as     pd
from   scipy.interpolate     import interp1d
from   scipy.ndimage.filters import convolve

class get_aux_data:

    def get_aux_data(self, path_mask:str=None, tell_lim:float=0.99, exclude_regions:list=None) -> None:
        """Get auxiliary data.

        :param path_mask: path to line mask (must be a CSV file where the wavelength column is "wave"), defaults to None
        :type path_mask: str, optional
        :param tell_lim: telluric depth limit, defaults to 0.99
        :type tell_lim: float, optional
        :param exclude_regions: wavelength intervals to be excluded, defaults to None
        :type exclude_regions: list, optional
        :return: None
        :rtype: None
        """
        
        # read wavelengh values
        wave_val = self.spec["wave_val"]

        # read resolution and medium
        resolution = self.spec["resolution"]
        medium     = self.spec["medium"]

        # read systematic velocity and maximum BERV
        vrad_sys = self.arve.star.stellar_parameters["vrad_sys"]
        berv_max = self.arve.star.stellar_parameters["berv_max"]

        # read constants
        c = self.arve.functions.constants["c"]

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

        # read self-provided mask
        if path_mask is not None:
            mask      = pd.read_csv(path_mask)
            mask_name = path_mask.split("/")[-1]
        # read closest mask
        else:
            mask      = pd.read_csv(path_aux_mask+file)
            mask_name = file
            if medium == "air":
                mask["wave"  ] = self.arve.functions.convert_vac_to_air(mask["wave"  ])
                mask["wave_l"] = self.arve.functions.convert_vac_to_air(mask["wave_l"])
                mask["wave_r"] = self.arve.functions.convert_vac_to_air(mask["wave_r"])
        mask_dict = {}
        for i in range(self.spec["Nord"]):
            idx_wave     = (mask.wave_l>=wave_val[i][0]) & (mask.wave_r<=wave_val[i][-1])
            mask_dict[i] = mask[idx_wave].reset_index(drop=True)

        # download closest spectrum from GitHub if it does not already exist in the package directory
        if not os.path.exists(path_aux_spec+file):
            path_github = "https://raw.githubusercontent.com/almoulla/arve/main/arve/aux_data/spectra/"+file
            urllib.request.urlretrieve(path_github, path_aux_spec)

        # read closest spectrum
        wave = pd.read_csv(path_aux_wave+"WAVE.csv.zip")["wave"]
        spec = pd.read_csv(path_aux_spec+file)
        if medium == "air":
            wave = self.arve.functions.convert_vac_to_air(wave)
        spec.insert(0, "wave", wave)

        # convolve spectrum
        if resolution is not None:
            spec["flux"] = _convolve_gaussian(spec.wave.to_numpy(), spec.flux.to_numpy(), resolution)

        # interpolate spectrum
        spec_wave = wave_val
        spec_flux = np.array([interp1d(spec.wave, spec.flux, kind="cubic", bounds_error=False)(spec_wave[i]) for i in range(self.spec["Nord"])])
        spec_temp = np.array([interp1d(spec.wave, spec.temp, kind="cubic", bounds_error=False)(spec_wave[i]) for i in range(self.spec["Nord"])])
        spec_dict = {
        "wave": spec_wave,
        "flux": spec_flux,
        "temp": spec_temp
        }

        # find left and right index of lines
        for i in range(self.spec["Nord"]):
            mask_dict[i].insert(mask_dict[i].columns.get_loc("wave_l")+2, "wave_l_idx", np.searchsorted(spec_dict["wave"][i], mask_dict[i].wave_l))
            mask_dict[i].insert(mask_dict[i].columns.get_loc("wave_r")+2, "wave_r_idx", np.searchsorted(spec_dict["wave"][i], mask_dict[i].wave_r))
            mask_dict[i]["wave_l"    ] = spec_dict["wave"][i][mask_dict[i].wave_l_idx]
            mask_dict[i]["wave_r"    ] = spec_dict["wave"][i][mask_dict[i].wave_r_idx]
        
        # read telluric spectrum
        wave = pd.read_csv(path_aux_wave+"WAVE.csv.zip")["wave"]
        tell = pd.read_csv(path_aux_tell+"TELL.csv.zip")
        if medium == "air":
            wave = self.arve.functions.convert_vac_to_air(wave)
        wave = self.arve.functions.doppler_shift(wave=wave, v=-vrad_sys)
        tell.insert(0, "wave", wave)

        # convolve telluric spectrum
        if resolution is not None:
            tell["flux"] = _convolve_gaussian(tell.wave.to_numpy(), tell.flux.to_numpy(), resolution)

        # interpolate telluric spectrum
        tell_wave = wave_val
        tell_flux = np.array([interp1d(tell.wave, tell.flux, kind="cubic", bounds_error=False)(tell_wave[i]) for i in range(self.spec["Nord"])])
        tell_dict = {
        "wave": tell_wave,
        "flux": tell_flux
        }

        # find left and right bounds of telluric bands
        tell_wave = np.array(tell["wave"])
        tell_flux = np.array(tell["flux"])
        idx_above = np.where(tell_flux>=tell_lim)[0]
        idx_below = np.where(tell_flux< tell_lim)[0]
        if len(idx_above) == 0:
            tell_wave_l = np.array([tell_wave[ 0]])
            tell_wave_r = np.array([tell_wave[-1]])
        elif len(idx_below) == 0:
            tell_wave_l = np.array([])
            tell_wave_r = np.array([])
        else:
            idx_l = idx_above[:-1][np.diff(idx_above)>1]
            idx_r = idx_below[:-1][np.diff(idx_below)>1]
            if tell_flux[ 0] <  tell_lim:
                idx_l = np.append(0, idx_l)
            if tell_flux[-1] >= tell_lim:
                idx_r = np.append(idx_r, idx_below[-1])
            if tell_flux[-1] <  tell_lim:
                idx_l = np.append(idx_l, idx_above[-1])
                idx_r = np.append(idx_r, len(tell_flux)-1)
            tell_wave_l = tell_wave[idx_l]*(1-berv_max/c)
            tell_wave_r = tell_wave[idx_r]*(1+berv_max/c)

        # create new telluric DataFrame with left and right bounds of telluric bands
        band = pd.DataFrame()
        band["wave_l"] = tell_wave_l
        band["wave_r"] = tell_wave_r

        # keep telluric bands which overlap with data
        band_dict = {}
        for i in range(self.spec["Nord"]):
            idx_wave     = (band.wave_l>=wave_val[i][-1]) | (band.wave_r<=wave_val[i][0]) == False
            band_dict[i] = band[idx_wave].reset_index(drop=True)
            band_dict[i]["wave_l"][band_dict[i]["wave_l"]<wave_val[i][ 0]] = wave_val[i][ 0]
            band_dict[i]["wave_r"][band_dict[i]["wave_r"]>wave_val[i][-1]] = wave_val[i][-1]

        # telluric criterion
        for i in range(self.spec["Nord"]):
            wave_l    = np.array(band_dict[i]["wave_l"])
            wave_r    = np.array(band_dict[i]["wave_r"])
            wave_band = np.array([wave_l,wave_r]).T
            mask_dict[i]["crit_tell"] = np.sum([(mask_dict[i]["wave"] > wave_band[j,0]) & (mask_dict[i]["wave"] < wave_band[j,1]) for j in range(len(wave_band))], axis=0).astype(bool) == False

        # exclude regions
        for i in range(self.spec["Nord"]):
            idx_excl = np.ones(len(mask_dict[i]), dtype=bool)
            if exclude_regions is not None:
                for j in range(len(mask_dict[i])):
                    if np.sum([(mask_dict[i]["wave"][j] > exclude_regions[k][0]) & (mask_dict[i]["wave"][j] < exclude_regions[k][1]) for k in range(len(exclude_regions))]) > 0:
                        idx_excl[j] = False
            mask_dict[i]["crit_excl"] = idx_excl

        # save
        self.aux_data = {"name": mask_name,
                         "mask": mask_dict,
                         "spec": spec_dict,
                         "tell": tell_dict,
                         "band": band_dict
                        }

        return None

def _convolve_gaussian(wave:list, flux:list, resolution:float) -> list:
    """Convolve spectrum with Gaussian instrumental profile. Adapted from: https://pysme-astro.readthedocs.io/en/latest/_modules/pysme/broadening.html

    :param wave: wavelength values
    :type wave: list
    :param flux: flux values
    :type flux: list
    :param resolution: instrumental resolution
    :type resolution: float
    :return: fluxed convolved with Gaussian instrumental profile
    :rtype: list
    """

    # nr. of wavelength points
    Nwave = len(wave)

    # half-width at half-maximum
    hwhm = 0.5*wave[0]/resolution

    # uniform dispersion
    wave_ran = wave[-1]-wave[0]
    wave_del = wave_ran/(Nwave-1)
    
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