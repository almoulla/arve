import numpy as     np
from   lmfit import Parameters, minimize

class fit_vpsd_coefficients:

    def fit_vpsd_coefficients(self) -> None:
        """Fit velocity power spectral density (VPSD) coefficients.

        :return: None
        :rtype: None
        """

        # read VPSD
        freq, vpsd, freq_avg, vpsd_avg = [self.vpsd[key] for key in ["freq", "vpsd", "freq_avg", "vpsd_avg"]]

        # LMFIT parameters
        params = Parameters()

        # loop components
        for comp in self.vpsd_components.keys():

            # component dictionary
            comp_dict = self.vpsd_components[comp]

            # coefficients and vary
            coef_val = comp_dict["coef_val"]
            vary     = comp_dict["vary"]

            # loop coefficients
            for i in range(len(coef_val)):

                # add parameters
                params.add(comp + "_" + str(i), value=coef_val[i], min=coef_val[i]/10, max=coef_val[i]*10, vary=vary[i])

        # fit coefficients
        c = minimize(_func_res, params, args=(self, freq_avg, vpsd_avg))

        # loop components
        for comp in self.vpsd_components.keys():

            # component dictionary
            comp_dict = self.vpsd_components[comp]

            # coefficients
            coef_val = comp_dict["coef_val"]
            coef_err = comp_dict["coef_err"]

            # loop coefficients
            for i in range(len(coef_val)):

                # update coefficients with fitted values
                coef_val[i] = c.params[comp + "_" + str(i)].value
                coef_err[i] = c.params[comp + "_" + str(i)].stderr
        
        return None

def _func_res(params, self, freq_avg, vpsd_avg):

    # empty array for sum of components
    vpsd_tot = np.zeros(len(freq_avg))

    # loop components
    for comp in self.vpsd_components.keys():

        # component dictionary
        comp_dict = self.vpsd_components[comp]

        # component type
        type = comp_dict["type"]

        # type Lorentz
        if type == "Lorentz":

            # unpack coefficients
            c0 = params[comp + "_0"]
            c1 = params[comp + "_1"]
            c2 = params[comp + "_2"]

            # compute component
            vpsd_comp = c0*c1**2/(c1**2+(freq_avg-c2)**2)

        # type Harvey
        if type == "Harvey":

            # unpack coefficients
            c0 = params[comp + "_0"]
            c1 = params[comp + "_1"]
            c2 = params[comp + "_2"]

            # compute component
            vpsd_comp = c0/(1+(c1*freq_avg)**c2)

        # type Constant
        if type == "Constant":
            
            # unpack coefficients
            c0 = params[comp + "_0"]

            # compute component
            vpsd_comp = c0

        # add component to sum
        vpsd_tot += vpsd_comp

    # logarithmic residual
    logres = np.log10(vpsd_avg) - np.log10(vpsd_tot)

    # return residual
    return logres