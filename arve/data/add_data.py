import glob
import pandas as pd

class add_data:

    def add_data(
        self,
        time_val      :list=None,
        vrad_val      :list=None,
        vrad_err      :list=None,
        wave_val      :list=None,
        flux_val      :list=None,
        flux_err      :list=None,
        resolution    :int =None,
        medium        :str="vac",
        path          :str= None,
        compression   :str= None,
        same_wave_grid:bool=False
        ) -> None:
        """Add data.

        :param time_val: time values, defaults to None
        :type time_val: list, optional
        :param vrad_val: radial velocity values, defaults to None
        :type vrad_val: list, optional
        :param vrad_err: radial velocity errors, defaults to None
        :type vrad_err: list, optional
        :param wave_val: wavelength values, defaults to None
        :type wave_val: list, optional
        :param flux_val: flux values, defaults to None
        :type flux_val: list, optional
        :param flux_err: flux errors, defaults to None
        :type flux_err: list, optional
        :param resolution: instrumental resolution, defaults to None
        :type resolution: int, optional
        :param medium: medium, defaults to "vac"
        :type medium: str, optional
        :param path: path to spectra which must be stored in CSV files, defaults to None
        :type path: str, optional
        :param compression: extension ('zip', 'gzip', etc) if the CSV files are compressed, defaults to None
        :type compression: str, optional
        :param same_wave_grid: all spectra already on the same wavelength grid, defaults to False
        :type same_wave_grid: bool, optional
        :return: None
        :rtype: None
        """
        
        # input from path
        if path is not None:
            
            # path to CSV files
            path_csv = path+'**/*.csv'

            # compressed extension
            if compression is not None:
                path_csv += '.'+compression
            
            # search files
            files = glob.glob(path_csv, recursive=True)
            files = sorted(files)

            # add first spectrum as reference
            df       = pd.read_csv(files[0])
            wave_val = df["wave_val"].to_numpy()
            flux_val = df["flux_val"].to_numpy()
            flux_val = flux_val.reshape(1,len(flux_val))
            flux_err = df["flux_err"].to_numpy()
            flux_err = flux_err.reshape(1,len(flux_err))
        
        else:

            # None entry for files
            files = None

        # add dictionaries with data
        self.time = {"time_val"      : time_val}
        self.vrad = {"vrad_val"      : vrad_val,
                     "vrad_err"      : vrad_err}
        self.spec = {"wave_val"      : wave_val,
                     "flux_val"      : flux_val,
                     "flux_err"      : flux_err,
                     "resolution"    : resolution,
                     "medium"        : medium,
                     "path"          : path,
                     "same_wave_grid": same_wave_grid,
                     "files"         : files}
        
        return None