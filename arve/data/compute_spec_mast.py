import numpy             as     np
import pandas            as     pd
from   scipy.interpolate import interp1d
import warnings
warnings.filterwarnings("ignore")

class compute_spec_mast:

    def compute_spec_mast(self, ofac=10) -> None:
        """Compute master spectrum.

        :param ofac: oversampling factor, defaults to 10
        :type ofac: int, optional
        :return: None
        :rtype: None
        """

        # read data from input
        if self.spec["path"] is None:

            # read data
            wave_val, flux_val, flux_err = [self.spec[key] for key in ["wave_val", "flux_val", "flux_err"]]

            # master spectrum
            mast_flux_val = np.average(flux_val, weights=1/flux_err**2, axis=0)
            
            # remove NaN
            mast_flux_val[np.isnan(mast_flux_val)] = 0

            # oversample
            func_flux     = interp1d(wave_val, mast_flux_val, "cubic")
            mast_wave_val = np.concatenate([np.linspace(wave_val[i], wave_val[i+1], ofac+1)[:ofac] for i in range(len(wave_val)-1)])
            mast_wave_val = np.append(mast_wave_val, wave_val[-1])
            mast_flux_val = func_flux(mast_wave_val)
        
        # read data from path
        else:

            # read data
            wave_val = self.spec["wave_val"]

            # nr. of spectra
            Nspec = len(self.spec["files"])

            # loop spectra
            for i in range(Nspec):

                # read CSV file
                df = pd.read_csv(self.spec["files"][i])
                
                # get flux and flux error if same wavelength grid
                if self.spec["same_wave_grid"]:
                    flux_val   = df["flux_val"].to_numpy()
                    flux_err   = df["flux_err"].to_numpy()
                
                # interpolate flux and flux error on reference wavelength grid
                else:
                    wave_val_i = df["wave_val"].to_numpy()
                    flux_val_i = df["flux_val"].to_numpy()
                    flux_err_i = df["flux_err"].to_numpy()
                    flux_val   = interp1d(wave_val_i, flux_val_i, kind="cubic", bounds_error=False)(wave_val)
                    flux_err   = interp1d(wave_val_i, flux_err_i, kind="cubic", bounds_error=False)(wave_val)
                
                # master spectrum
                if i == 0:
                    mast_flux_val = flux_val
                    mast_flux_err = flux_err
                else:
                    flux_val_arr = np.array([mast_flux_val,flux_val])
                    flux_err_arr = np.array([mast_flux_err,flux_err])
                    mast_flux_val = np.average(flux_val_arr, weights=1/flux_err_arr**2, axis=0)
                    mast_flux_err = np.sqrt(1/np.sum(1/flux_err_arr**2, axis=0))

            # remove NaN
            mast_flux_val[np.isnan(mast_flux_val)] = 0
            
            # oversample
            func_flux     = interp1d(wave_val, mast_flux_val, "cubic")
            mast_wave_val = np.concatenate([np.linspace(wave_val[i], wave_val[i+1], ofac+1)[:ofac] for i in range(len(wave_val)-1)])
            mast_wave_val = np.append(mast_wave_val, wave_val[-1])
            mast_flux_val = func_flux(mast_wave_val)

        # save spectral data
        self.spec_mast = {"wave_val": mast_wave_val, "flux_val": mast_flux_val}

        return None