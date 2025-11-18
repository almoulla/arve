import numpy             as     np
from   scipy.interpolate import interp1d
from   scipy.signal      import savgol_filter
import warnings
warnings.filterwarnings("ignore")

class compute_spec_reference:

    def compute_spec_reference(
        self,
        N_spec   : int | None = None ,
        oversamp : int        = 10   ,
        smooth   : bool       = False,
        ) -> None:
        """Compute reference spectrum.

        Parameters
        ----------
        N_spec : int | None, optional
            number of spectra used to build the reference (if None, all input spectra will be considered), by default None
        oversamp : int, optional
            oversampling factor of the wavelength grid, by default 10
        smooth : bool, optional
            smooth the individual spectra before building the reference, by default False

        Returns
        -------
        None
            None
        """

        # read data
        wave_val = self.spec["wave_val"]
        N_ord    = self.spec["N_ord"]
        N_pix    = self.spec["N_pix"]
        kind     = self.spec["interpolation"]

        # nr. of spectra
        if N_spec is None:
            N_spec = self.spec["N_spec"]

        # read data from input
        if self.spec["path"] is None:

            # read data
            flux_val, flux_err = [self.spec[key] for key in ["flux_val", "flux_err"]]

            # reference spectrum
            ref_wave_val = wave_val
            ref_flux_val = np.array([np.average(flux_val[:N_spec,i], weights=1/flux_err[:N_spec,i]**2, axis=0) for i in range(N_ord)])
            ref_flux_err = np.array([np.sqrt(1/np.sum(1/flux_err[:N_spec,i]**2, axis=0))                       for i in range(N_ord)])
        
        # read data from path
        else:

            # loop spectra
            for i in range(N_spec):

                # read spectrum
                _, flux_val, flux_err = self.read_spec(i)

                # smooth spectrum
                if smooth:
                    for j in range(N_ord):
                        flux_val[j] = savgol_filter(flux_val[j], window_length=15, polyorder=2)
                
                # reference spectrum
                if i == 0:
                    ref_wave_val = wave_val
                    ref_flux_val = flux_val
                    ref_flux_err = flux_err
                else:
                    flux_val_arr = np.array([ref_flux_val,flux_val])
                    flux_err_arr = np.array([ref_flux_err,flux_err])
                    ref_flux_val = np.array([np.average(flux_val_arr[:,j], weights=1/flux_err_arr[:,j]**2, axis=0) for j in range(N_ord)])
                    ref_flux_err = np.array([np.sqrt(1/np.sum(1/flux_err_arr[:,j]**2, axis=0))                     for j in range(N_ord)])

        # oversample
        ref_wave_val_samp = np.zeros((N_ord,N_pix+(oversamp-1)*(N_pix-1)))
        ref_flux_val_samp = np.zeros((N_ord,N_pix+(oversamp-1)*(N_pix-1)))
        ref_flux_err_samp = np.zeros((N_ord,N_pix+(oversamp-1)*(N_pix-1)))
        for i in range(N_ord):
            ref_wave_val_samp[i] = np.append(np.concatenate([np.linspace(ref_wave_val[i][j], ref_wave_val[i][j+1], oversamp+1)[:oversamp] for j in range(N_pix-1)]), ref_wave_val[i][-1])
            idx_nan_samp = np.isnan(ref_wave_val_samp[i])
            idx_val_samp = ~idx_nan_samp
            if np.sum(idx_val_samp) > 0:
                idx_val = ~np.isnan(ref_wave_val[i])
                ref_flux_val_samp[i,idx_val_samp] = interp1d(ref_wave_val[i,idx_val], ref_flux_val[i,idx_val], kind=kind)(ref_wave_val_samp[i,idx_val_samp])
                ref_flux_err_samp[i,idx_val_samp] = interp1d(ref_wave_val[i,idx_val], ref_flux_err[i,idx_val], kind=kind)(ref_wave_val_samp[i,idx_val_samp])
            if np.sum(idx_nan_samp) > 0:
                ref_flux_val_samp[i,idx_nan_samp] = np.nan
                ref_flux_err_samp[i,idx_nan_samp] = np.nan

        # save spectral data
        self.spec_reference = {"wave_val": ref_wave_val_samp,
                               "flux_val": ref_flux_val_samp,
                               "flux_err": ref_flux_err_samp}

        return None