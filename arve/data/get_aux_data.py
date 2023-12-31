import numpy                 as     np
import os
import pandas                as     pd
import pkg_resources
from   scipy                 import ndimage
from   scipy.interpolate     import interp1d
from   scipy.ndimage.filters import convolve
from   scipy.signal          import argrelextrema

class get_aux_data:

    def get_aux_data(self, path_mask:str=None, tell_lim:float=0.99) -> None:
        """Get auxiliary data.

        :param path_mask: path to line mask (must be a CSV file where the wavelength column is "wave"), defaults to None
        :type path_mask: str, optional
        :param tell_lim: telluric depth limit, defaults to 0.99
        :type tell_lim: float, optional
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

        # path to auxiliary data
        path_aux_data = pkg_resources.resource_filename("arve", "aux_data/")

        # search masks
        masks            = sorted(os.listdir(path_aux_data+"masks/"))
        masks            = [mask for mask in masks if mask.endswith(".csv.zip")]
        sptype_masks     = [mask.split(".")[0] for mask in masks]
        sptype_num_masks = np.array([self.arve.functions.sptype_to_num(sptype=sptype) for sptype in sptype_masks])

        # search spec
        specs            = sorted(os.listdir(path_aux_data+"spectra/"))
        specs            = [spec for spec in specs if spec.endswith(".csv.zip")]
        sptype_specs     = [spec.split(".")[0] for spec in specs]
        sptype_num_specs = np.array([self.arve.functions.sptype_to_num(sptype=sptype) for sptype in sptype_specs])

        # spectral type as number
        sptype_num       = self.arve.functions.sptype_to_num(sptype=self.arve.star.stellar_parameters["sptype"])

        # read self-provided mask
        if path_mask is not None:
            mask      = pd.read_csv(path_mask)
            mask_name = path_mask.split("/")[-1]
        # read closest mask
        else:
            idx_mask  = np.argmin(np.abs(sptype_num-sptype_num_masks))
            mask      = pd.read_csv(path_aux_data+"masks/"+masks[idx_mask])
            mask_name = masks[idx_mask]
            if medium == "air":
                mask["wave"  ] = self.arve.functions.convert_vac_to_air(mask["wave"  ])
                mask["wave_l"] = self.arve.functions.convert_vac_to_air(mask["wave_l"])
                mask["wave_r"] = self.arve.functions.convert_vac_to_air(mask["wave_r"])
        mask_dict = {}
        for i in range(self.spec["Nord"]):
            idx_wave     = (mask.wave_l>=wave_val[i][0]) & (mask.wave_r<=wave_val[i][-1])
            mask_dict[i] = mask[idx_wave].reset_index(drop=True)

        # read closest spectrum
        idx_spec = np.argmin(np.abs(sptype_num-sptype_num_specs))
        wave     = pd.read_csv(path_aux_data+"wavelengths/WAVE.csv.zip")["wave"]
        spec     = pd.read_csv(path_aux_data+"spectra/"+specs[idx_spec])
        if medium == "air":
            wave = self.arve.functions.convert_vac_to_air(wave)
        spec.insert(0, "wave", wave)

        # convolve spectrum
        if resolution is not None:
            spec["flux"] = _convolve_gaussian(spec.wave.to_numpy(), spec.flux.to_numpy(), resolution)

        # interpolate spectrum
        spec_wave = wave_val
        spec_flux = np.array([interp1d(spec.wave, spec.flux, kind="cubic")(spec_wave[i]) for i in range(self.spec["Nord"])])
        spec_temp = np.array([interp1d(spec.wave, spec.temp, kind="cubic")(spec_wave[i]) for i in range(self.spec["Nord"])])
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
        wave = pd.read_csv(path_aux_data+"wavelengths/WAVE.csv.zip")["wave"]
        tell = pd.read_csv(path_aux_data+"tellurics/TELL.csv.zip")
        if medium == "air":
            wave = self.arve.functions.convert_vac_to_air(wave)
        wave = self.arve.functions.doppler_shift(wave=wave, v=-vrad_sys)
        tell.insert(0, "wave", wave)

        # convolve telluric spectrum
        if resolution is not None:
            tell["flux"] = _convolve_gaussian(tell.wave.to_numpy(), tell.flux.to_numpy(), resolution)

        # interpolate telluric spectrum
        tell_wave = wave_val
        tell_flux = np.array([interp1d(tell.wave, tell.flux, kind="cubic")(tell_wave[i]) for i in range(self.spec["Nord"])])
        tell_dict = {
        "wave": tell_wave,
        "flux": tell_flux
        }

        # find telluric lines
        tell_wave = np.array(tell["wave"])
        tell_flux = np.array(tell["flux"])
        idx = argrelextrema(tell_flux, np.less)[0]
        wc  = tell_wave[idx]
        fc  = tell_flux[idx]
        idx = fc < tell_lim
        wc  = wc[idx]
        fc  = fc[idx]

        # wavelength window around telluric lines due to BERV
        dwave = wc*berv_max/c

        # left and right bounds of telluric lines
        wl = wc-dwave
        wr = wc+dwave

        # label overlapping telluric lines
        label_map, Nlabel = ndimage.label(wl[1:] < wr[:-1])

        # empty arrays for left and right bounds of telluric bands
        wave_l = np.zeros(Nlabel)
        wave_r = np.zeros(Nlabel)

        # loop labels
        for i in range(Nlabel):

            # find left- and rightmost edges of overlapping telluric lines
            idx = np.where(label_map == i+1)[0]
            wave_l[i] = wl[idx[ 0]  ]
            wave_r[i] = wr[idx[-1]+1]
            
            # label last telluric line in each band
            if idx[-1]+1 < len(label_map):
                label_map[idx[-1]+1] = i+1

        # nr. of single telluric lines
        Nsingle = np.sum(label_map == 0)

        # loop single telluric lines
        for i in range(Nsingle):

            # find index of single telluric line
            idx = np.where(label_map == 0)[0][i]

            # append telluric line bounds to list of telluric bands
            wave_l = np.append(wave_l, wl[idx])
            wave_r = np.append(wave_r, wr[idx])

        # sort telluric bounds
        idx = np.argsort(wave_l)
        wave_l = wave_l[idx]
        wave_r = wave_r[idx]

        # create new telluric DataFrame with left and right edges of telluric bands
        band = pd.DataFrame()
        band["wave_l"] = wave_l
        band["wave_r"] = wave_r

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

        # save
        self.aux_data = {"name": mask_name,
                         "mask": mask_dict,
                         "spec": spec_dict,
                         "tell": tell_dict,
                         "band": band_dict
                        }

        return None

def _convolve_gaussian(wave:list, flux:list, resolution:float) -> list:
    """Convolve spectrum with Gaussian instrumental profile.

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