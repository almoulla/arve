from   astropy.io        import fits
import numpy             as     np
import pandas            as     pd
from   scipy.interpolate import interp1d

class read_spec:

    def read_spec(self, i:int) -> tuple:
        """Read spectrum.

        :param i: spectrum index
        :type i: int
        :return: wavelength values, flux values and flux errors of i:th spectrum
        :rtype: tuple
        """

        # read data
        vrad_sys = self.arve.star.stellar_parameters["vrad_sys"]

        # read data from input
        if self.spec["path"] is None:
            wave_val = self.spec["wave_val"]
            flux_val = self.spec["flux_val"][i]
            flux_err = self.spec["flux_err"][i]
        
        # read data from path
        if self.spec["path"] is not None:

            # read file: CSV
            if self.spec["extension"] == "csv":
                file = pd.read_csv(self.spec["files"][i])
            
            # read file: NPZ
            if self.spec["extension"] == "npz":
                file = np.load(self.spec["files"][i])
                self.time["time_val"][i] = file["time_val"]

            # read file: FITS
            if self.spec["extension"] == "fits":

                # instrument: NIRPS
                if self.spec["instrument"] == "nirps":
                    if i == 0: self.spec["medium"] = "air"
                    hdul = fits.open(self.spec["files"][i])
                    self.time["time_val"][i] = hdul[0].header["HIERARCH ESO QC BJD"]
                    file = {"wave_val": hdul[5].data,
                            "flux_val": hdul[1].data,
                            "flux_err": hdul[2].data}
                    hdul.close()
            
            # shift wavelengths and reshape arrays
            wave_val = np.array(file["wave_val"])
            flux_val = np.array(file["flux_val"])
            flux_err = np.array(file["flux_err"])
            wave_val = self.arve.functions.doppler_shift(wave=wave_val, v=-vrad_sys)
            if self.spec["format"] == "s1d":
                wave_val = wave_val.reshape(1,wave_val.shape[0])
                flux_val = flux_val.reshape(1,flux_val.shape[0])
                flux_err = flux_err.reshape(1,flux_err.shape[0])
            
            # interpolate flux values and errors on reference wavelength grid
            if (self.spec["same_wave_grid"] == False) & (i > 0):
                wave_val_inter = self.spec["wave_val"]
                flux_val_inter = np.zeros_like(wave_val_inter)
                flux_err_inter = np.zeros_like(wave_val_inter)
                for j in range(self.spec["Nord"]):
                    flux_val_inter[j] = interp1d(wave_val[j], flux_val[j], kind="cubic", bounds_error=False)(wave_val_inter[j])
                    flux_err_inter[j] = interp1d(wave_val[j], flux_err[j], kind="cubic", bounds_error=False)(wave_val_inter[j])
                wave_val = wave_val_inter
                flux_val = flux_val_inter
                flux_err = flux_err_inter
        
        return wave_val, flux_val, flux_err