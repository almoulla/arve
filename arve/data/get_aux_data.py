import numpy                 as     np
import os
import pandas                as     pd
import pkg_resources
from   scipy.interpolate     import interp1d
from   scipy.ndimage.filters import convolve

class get_aux_data:

    def get_aux_data(self) -> None:
        """Get auxiliary data.

        :return: None
        :rtype: None
        """
        
        # read wavelengh range
        wave_val = self.spec["wave_val"]
        wave_min = min(wave_val)
        wave_max = max(wave_val)

        # read resolution and medium
        resolution = self.spec["resolution"]
        medium     = self.spec["medium"]

        # path to auxiliary data
        path_aux_data = pkg_resources.resource_filename("arve", "aux_data/")

        # search masks
        masks            = os.listdir(path_aux_data+"masks/")
        masks            = [mask for mask in masks if mask.endswith(".csv.zip")]
        sptype_masks     = [mask.split(".")[0] for mask in masks]
        sptype_num_masks = np.array([self.arve.functions.sptype_to_num(sptype=sptype) for sptype in sptype_masks])

        # search spec
        specs            = os.listdir(path_aux_data+"spectra/")
        specs            = [spec for spec in specs if spec.endswith(".csv.zip")]
        sptype_specs     = [spec.split(".")[0] for spec in specs]
        sptype_specs.remove('WAVE')
        sptype_num_specs = np.array([self.arve.functions.sptype_to_num(sptype=sptype) for sptype in sptype_specs])

        # spectral type as number
        sptype_num       = self.arve.functions.sptype_to_num(sptype=self.arve.star.stellar_parameters["sptype"])

        # read closest mask
        idx_mask = np.argmin(np.abs(sptype_num-sptype_num_masks))
        mask     = pd.read_csv(path_aux_data+"masks/"+masks[idx_mask])
        if medium == "air":
            mask["wave"  ] = self.arve.functions.convert_vac_to_air(mask["wave"  ])
            mask["wave_l"] = self.arve.functions.convert_vac_to_air(mask["wave_l"])
            mask["wave_r"] = self.arve.functions.convert_vac_to_air(mask["wave_r"])
        idx_wave = (mask.wave_l>=wave_min) & (mask.wave_r<=wave_max)
        mask     = mask[idx_wave].reset_index(drop=True)

        # read closest spec
        idx_spec = np.argmin(np.abs(sptype_num-sptype_num_specs))
        spec     = pd.read_csv(path_aux_data+"spectra/"+specs[idx_spec])
        wave     = pd.read_csv(path_aux_data+"spectra/WAVE.csv.zip")["wave"]
        if medium == "air":
            wave = self.arve.functions.convert_vac_to_air(wave)
        spec.insert(0, "wave", wave)

        # interpolate spectrum
        spec_new         = pd.DataFrame()
        spec_new["wave"] = wave_val
        spec_new["flux"] = interp1d(spec.wave, spec.flux, kind="cubic")(spec_new.wave)
        spec_new["temp"] = interp1d(spec.wave, spec.temp, kind="cubic")(spec_new.wave)
        spec             = spec_new

        # convolve spectrum
        if resolution is not None:
            spec["flux"] = _convolve_gaussian(spec.wave.to_numpy(), spec.flux.to_numpy(), resolution)

        # read tellurics
        tell     = pd.read_csv(path_aux_data+"tellurics/TELL.csv.zip")
        tell["wave_l"] = self.arve.functions.doppler_shift(tell["wave_l"], -self.arve.star.stellar_parameters["vrad_sys"])
        tell["wave_r"] = self.arve.functions.doppler_shift(tell["wave_r"], -self.arve.star.stellar_parameters["vrad_sys"])
        if medium == "air":
            tell["wave_l"] = self.arve.functions.convert_vac_to_air(tell["wave_l"])
            tell["wave_r"] = self.arve.functions.convert_vac_to_air(tell["wave_r"])
        idx_wave = (tell.wave_l>=wave_min) & (tell.wave_r<=wave_max)
        tell     = tell[idx_wave].reset_index(drop=True)

        # find left and right index of lines
        mask.insert(mask.columns.get_loc("wave_l")+2, "wave_l_idx", np.searchsorted(spec.wave, mask.wave_l))
        mask.insert(mask.columns.get_loc("wave_r")+2, "wave_r_idx", np.searchsorted(spec.wave, mask.wave_r))
        mask["wave_l"    ] = spec.wave[mask.wave_l_idx].to_numpy()
        mask["wave_r"    ] = spec.wave[mask.wave_r_idx].to_numpy()

        # save
        self.aux_data = {"name": masks[idx_mask],
                         "mask": mask,
                         "spec": spec,
                         "tell": tell}

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