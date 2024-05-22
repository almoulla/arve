import numpy             as     np
from   scipy.interpolate import interp1d
import warnings
warnings.filterwarnings("ignore")

class compute_spec_mast:

    def compute_spec_mast(self, Nspec:int=None, ofac:int=10) -> None:
        """Compute master spectrum.

        :param Nspec: number of spectra used to build the master (if None, all input spectra will be considered), defaults to None
        :type Nspec: int, optional
        :param ofac: oversampling factor, defaults to 10
        :type ofac: int, optional
        :return: None
        :rtype: None
        """

        # read data
        wave_val = self.spec["wave_val"]
        Nord     = self.spec["Nord"]
        Npix     = self.spec["Npix"]

        # read data from input
        if self.spec["path"] is None:

            # read data
            flux_val, flux_err = [self.spec[key] for key in ["flux_val", "flux_err"]]

            # master spectrum
            mast_wave_val = wave_val
            mast_flux_val = np.array([np.average(flux_val[:,i], weights=1/flux_err[:,i]**2, axis=0) for i in range(Nord)])
            mast_flux_err = np.array([np.sqrt(1/np.sum(1/flux_err[:,i]**2, axis=0))                 for i in range(Nord)])
        
        # read data from path
        else:

            # nr. of spectra
            if Nspec is None:
                Nspec = len(self.spec["files"])

            # loop spectra
            for i in range(Nspec):

                # read spectrum
                _, flux_val, flux_err = self.read_spec(i)
                
                # master spectrum
                if i == 0:
                    mast_wave_val = wave_val
                    mast_flux_val = flux_val
                    mast_flux_err = flux_err
                else:
                    flux_val_arr  = np.array([mast_flux_val,flux_val])
                    flux_err_arr  = np.array([mast_flux_err,flux_err])
                    mast_flux_val = np.array([np.average(flux_val_arr[:,j], weights=1/flux_err_arr[:,j]**2, axis=0) for j in range(Nord)])
                    mast_flux_err = np.array([np.sqrt(1/np.sum(1/flux_err_arr[:,j]**2, axis=0))                     for j in range(Nord)])

        # remove NaN
        idx = np.isnan(mast_flux_val) | np.isnan(mast_flux_err)
        mast_flux_val[idx] = 0
        mast_flux_err[idx] = 0

        # oversample
        mast_wave_val_samp = np.zeros((Nord,Npix+(ofac-1)*(Npix-1)))
        mast_flux_val_samp = np.zeros((Nord,Npix+(ofac-1)*(Npix-1)))
        mast_flux_err_samp = np.zeros((Nord,Npix+(ofac-1)*(Npix-1)))
        for i in range(Nord):
            func_flux_val         = interp1d(mast_wave_val[i], mast_flux_val[i], "cubic")
            func_flux_err         = interp1d(mast_wave_val[i], mast_flux_err[i], "cubic")
            mast_wave_val_samp[i] = np.append(np.concatenate([np.linspace(mast_wave_val[i][j], mast_wave_val[i][j+1], ofac+1)[:ofac] for j in range(Npix-1)]), mast_wave_val[i][-1])
            mast_flux_val_samp[i] = func_flux_val(mast_wave_val_samp[i])
            mast_flux_err_samp[i] = func_flux_err(mast_wave_val_samp[i])

        # save spectral data
        self.spec_mast = {"wave_val": mast_wave_val_samp,
                          "flux_val": mast_flux_val_samp,
                          "flux_err": mast_flux_err_samp}

        return None