from   lmfit import Parameters, minimize
import numpy as     np

class fit_vpsd_coefficients:

    def fit_vpsd_coefficients(
        self,
        coef_bound : float = 10,
        ) -> None:
        """Fit velocity power spectral density (VPSD) coefficients.

        Parameters
        ----------
        coef_bound : float, optional
            multiplicative factor within which the fitted coefficient are bounded with respect to their theoretical values, by default 10

        Returns
        -------
        None
            None
        """

        # read VPSD
        freq, vpsd, freq_avg, vpsd_avg = [self.vpsd[key] for key in ["freq", "vpsd", "freq_avg", "vpsd_avg"]]

        # LMFIT parameters
        params = Parameters()

        # loop components
        for comp_name in self.vpsd_components.keys():

            # component dictionary
            comp_dict = self.vpsd_components[comp_name]

            # coefficients and vary
            coef_val  = comp_dict["coef_val"]
            coef_vary = comp_dict["coef_vary"]

            # loop coefficients
            for i in range(len(coef_val)):

                # add parameters
                params.add(comp_name + "_" + str(i), value=coef_val[i], min=coef_val[i]/coef_bound, max=coef_val[i]*coef_bound, vary=coef_vary[i])

        # fit coefficients
        coef_fit = minimize(_func_res, params, args=(self, freq_avg, vpsd_avg))

        # loop components
        for comp_name in self.vpsd_components.keys():

            # component dictionary
            comp_dict = self.vpsd_components[comp_name]

            # coefficients
            coef_val = comp_dict["coef_val"]
            coef_err = comp_dict["coef_err"]

            # loop coefficients
            for i in range(len(coef_val)):

                # update coefficients with fitted values
                coef_val[i] = coef_fit.params[comp_name + "_" + str(i)].value
                coef_err[i] = coef_fit.params[comp_name + "_" + str(i)].stderr
        
        return None

def _func_res(
    params   : dict[str,float],
    self,
    freq_avg : np.ndarray     ,
    vpsd_avg : np.ndarray
    ) -> np.ndarray:
    """VPSD residual function.

    Parameters
    ----------
    params : dict[str,float]
        VPSD coefficients
    freq_avg : np.ndarray
        Binned frequency
    vpsd_avg : np.ndarray
        Binned VPSD

    Returns
    -------
    np.ndarray
        Logarithmic residuals between binned VPSD and modeled VPSD (sum of all included components).
    """

    # empty array for sum of components
    vpsd_tot = np.zeros(len(freq_avg))

    # loop components
    for comp_name in self.vpsd_components.keys():

        # component dictionary
        comp_dict = self.vpsd_components[comp_name]

        # function type
        func_type = comp_dict["func_type"]

        # type Lorentz
        if func_type == "lorentz":

            # unpack coefficients
            c0 = params[comp_name + "_0"]
            c1 = params[comp_name + "_1"]
            c2 = params[comp_name + "_2"]

            # compute component
            vpsd_comp = c0*c1**2/(c1**2+(freq_avg-c2)**2)

        # type Harvey
        if func_type == "harvey":

            # unpack coefficients
            c0 = params[comp_name + "_0"]
            c1 = params[comp_name + "_1"]
            c2 = params[comp_name + "_2"]

            # compute component
            vpsd_comp = c0/(1+(c1*freq_avg)**c2)

        # type constant
        if func_type == "constant":
            
            # unpack coefficients
            c0 = params[comp_name + "_0"]

            # compute component
            vpsd_comp = c0

        # add component to sum
        vpsd_tot += vpsd_comp

    # logarithmic residual
    logres = np.log10(vpsd_avg) - np.log10(vpsd_tot)

    # return residual
    return logres