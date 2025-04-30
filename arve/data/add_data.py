import glob
import numpy as np

from typing import Literal

class add_data:

    def add_data(
        self,
        time_val       : np.ndarray                    | None = None ,
        vrad_val       : np.ndarray                    | None = None ,
        vrad_err       : np.ndarray                    | None = None ,
        berv_val       : np.ndarray                    | None = None ,
        wave_val       : np.ndarray                    | None = None ,
        flux_val       : np.ndarray                    | None = None ,
        flux_err       : np.ndarray                    | None = None ,
        medium         : Literal["vac", "air"]                = "vac",
        format         : Literal["s1d", "s2d"]                = "s1d",
        path           : str                           | None = None ,
        extension      : Literal["fits", "npz", "csv"] | None = None ,
        compression    : str                           | None = None ,
        instrument     : Literal["espresso",
                                 "harps",
                                 "harps-n",
                                 "kpf",
                                 "neid",
                                 "nirps",
                                 "spirou"]             | None = None ,
        resolution     : float                         | None = None ,
        berv_corrected : bool                                 = True ,
        same_wave_grid : bool                                 = False
        ) -> None:
        """Add data.

        Parameters
        ----------
        time_val : np.ndarray | None, optional
            time values, by default None
        vrad_val : np.ndarray | None, optional
            radial velocity values, by default None
        vrad_err : np.ndarray | None, optional
            radial velocity errors, by default None
        berv_val : np.ndarray | None, optional
            barycentric-Earth radial velocity (BERV) values, by default None
        wave_val : np.ndarray | None, optional
            wavelength values, by default None
        flux_val : np.ndarray | None, optional
            flux values, by default None
        flux_err : np.ndarray | None, optional
            flux errors, by default None
        medium : Literal[&quot;vac&quot;, &quot;air&quot;], optional
            medium of recorded wavelengths, by default "vac"
        format : Literal[&quot;s1d&quot;, &quot;s2d&quot;], optional
            spectral format, by default "s1d"
        path : str | None, optional
            path to spectra, by default None
        extension : Literal[&quot;fits&quot;, &quot;npz&quot;, &quot;csv&quot;] | None, optional
            extension of files containing the spectra, by default None
        compression : str | None, optional
            compression of files containing the spectra, by default None
        instrument : Literal[&quot;espresso&quot;, &quot;harps&quot;, &quot;harps-n&quot;, &quot;kpf&quot;, &quot;neid&quot;, &quot;nirps&quot;, &quot;spirou&quot;] | None, optional
            instrument name (if the spectra are in FITS files), by default None
        resolution : float | None, optional
            instrumental resolution (if the spectra are in NPZ or CSV files), by default None
        berv_corrected : bool, optional
            spectra already BERV-corrected, by default True
        same_wave_grid : bool, optional
            spectra already on the same wavelength grid, by default False

        Returns
        -------
        None
            None
        """
        
        # read data
        if self.arve.star.stellar_parameters is not None:
            vrad_sys = self.arve.star.stellar_parameters["vrad_sys"]

        # input from keyword arguments
        if path is None:

            # shift wavelength grid to system restframe
            if wave_val is not None:
                wave_val = self.arve.functions.doppler_shift(wave=wave_val, v=-vrad_sys)

            # None entry for files
            files = None

        # input from path
        if path is not None:
            
            # path to files
            path_ext = path+f"**/*.{extension}"
            
            # path with compression
            if compression is not None:
                path_ext += "."+compression
            
            # search files
            files = glob.glob(path_ext, recursive=True)
            files = sorted(files)

        # add dictionaries with data
        self.time = {
            "time_val"      : time_val,
            "berv_val"      : berv_val
            }
        self.vrad = {
            "vrad_val"      : vrad_val,
            "vrad_err"      : vrad_err
            }
        self.spec = {
            "wave_val"      : wave_val,
            "flux_val"      : flux_val,
            "flux_err"      : flux_err,
            "format"        : format,
            "resolution"    : resolution,
            "medium"        : medium,
            "path"          : path,
            "extension"     : extension,
            "instrument"    : instrument,
            "berv_corrected": berv_corrected,
            "same_wave_grid": same_wave_grid,
            "files"         : files
            }

        # if input from path, make special additions
        if path is not None:

            # if not provided, make time values an array with zeros to be populated
            if time_val is None:
                self.time["time_val"] = np.zeros(len(files))
            # if not provided, make BERV values an array with zeros to be populated
            if berv_val is None:
                self.time["berv_val"] = np.zeros(len(files))
            
            # read reference (0th) spectrum
            wave_val, flux_val, flux_err = self.read_spec(0)

            # reshape arrays
            flux_val = flux_val.reshape(1,flux_val.shape[0],flux_val.shape[1])
            flux_err = flux_err.reshape(1,flux_err.shape[0],flux_err.shape[1])

            # add reference spectrum
            self.spec["wave_val"] = wave_val
            self.spec["flux_val"] = flux_val
            self.spec["flux_err"] = flux_err

        # add nr. of spectra, orders and pixels
        if wave_val is not None:
            if path is None:
                if flux_val is not None:
                    self.spec["N_spec"] = flux_val.shape[0]
            else:
                self.spec["N_spec"] = len(files)
            self.spec["N_ord"] = wave_val.shape[0]
            self.spec["N_pix"] = wave_val.shape[1]
        
        return None