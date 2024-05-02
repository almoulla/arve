import glob
import numpy as np

class add_data:

    def add_data(
        self,
        time_val      :list=None ,
        vrad_val      :list=None ,
        vrad_err      :list=None ,
        berv_val      :list=None ,
        wave_val      :list=None ,
        flux_val      :list=None ,
        flux_err      :list=None ,
        format        :str="s1d" ,
        resolution    :int =None ,
        medium        :str="vac" ,
        path          :str= None ,
        extension     :str="csv" ,
        compression   :str= None ,
        instrument    :str= None ,
        berv_corrected:bool=True ,
        same_wave_grid:bool=False,
        ) -> None:
        """Add data.

        :param time_val: time values, defaults to None
        :type time_val: list, optional
        :param vrad_val: radial velocity values, defaults to None
        :type vrad_val: list, optional
        :param vrad_err: radial velocity errors, defaults to None
        :type vrad_err: list, optional
        :param berv_val: BERV values, defaults to None
        :type berv_val: list, optional
        :param wave_val: wavelength values, defaults to None
        :type wave_val: list, optional
        :param flux_val: flux values, defaults to None
        :type flux_val: list, optional
        :param flux_err: flux errors, defaults to None
        :type flux_err: list, optional
        :param format: spectral format (S1D or S2D), defaults to "s1d"
        :type format: str, optional
        :param resolution: instrumental resolution, defaults to None
        :type resolution: int, optional
        :param medium: medium, defaults to "vac"
        :type medium: str, optional
        :param path: path to spectra which must be stored in CSV, NPZ or FITS files, defaults to None
        :type path: str, optional
        :param extension: extension (csv, npz or fits), defaults to "csv"
        :type extension: str, optional
        :param compression: compression (zip, gzip, etc) if the CSV files are compressed, defaults to None
        :type compression: str, optional
        :param instrument: instrument name if the spectra are in FITS files, defaults to None
        :type instrument: str, optional
        :param berv_corrected: spectra already BERV corrected, defaults to True
        :type berv_corrected: bool, optional
        :param same_wave_grid: spectra already on the same wavelength grid, defaults to False
        :type same_wave_grid: bool, optional
        :return: None
        :rtype: None
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
                    self.spec["Nspec"] = flux_val.shape[0]
            else:
                self.spec["Nspec"] = len(files)
            self.spec["Nord"] = wave_val.shape[0]
            self.spec["Npix"] = wave_val.shape[1]
        
        return None