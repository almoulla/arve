import numpy as np
from lmfit import Parameters, minimize

from star import Star


def fit_vpsd_coefficients(self: Star) -> None:
    """Fit velocity power spectral density (VPSD) coefficients.

    :return: None
    """
    # read VPSD
    freq, vpsd, freq_avg, vpsd_avg = (
        self.vpsd[key] for key in ["freq", "vpsd", "freq_avg", "vpsd_avg"]
    )

    # LMFIT parameters
    params = Parameters()

    # loop components
    for comp in self.vpsd_components:
        # component dictionary
        comp_dict = self.vpsd_components[comp]

        # coefficients and vary
        coef_val = comp_dict["coef_val"]
        vary = comp_dict["vary"]

        # loop coefficients
        for i in range(len(coef_val)):
            # add parameters
            params.add(
                f"{comp}_{str(i)}",
                value=coef_val[i],
                min=coef_val[i] / 10,
                max=coef_val[i] * 10,
                vary=vary[i],
            )

    # fit coefficients
    c = minimize(_func_res, params, args=(self, freq_avg, vpsd_avg))

    # loop components
    for comp in self.vpsd_components:
        # component dictionary
        comp_dict = self.vpsd_components[comp]

        # coefficients
        coef_val = comp_dict["coef_val"]
        coef_err = comp_dict["coef_err"]

        # loop coefficients
        for i in range(len(coef_val)):
            # update coefficients with fitted values
            coef_val[i] = c.params[f"{comp}_{str(i)}"].value
            coef_err[i] = c.params[f"{comp}_{str(i)}"].stderr


def _func_res(
    params: Parameters, self: Star, freq_avg: np.ndarray, vpsd_avg: np.ndarray
) -> np.ndarray:
    # empty array for sum of components
    vpsd_tot = np.zeros(len(freq_avg))

    # loop components
    for comp in self.vpsd_components:
        # component dictionary
        comp_dict = self.vpsd_components[comp]

        # component type
        component_type = comp_dict["type"]

        # unpack coefficients
        c0 = params[f"{comp}_0"]
        c1 = params[f"{comp}_1"]
        c2 = params[f"{comp}_2"]

        # type Constant
        if component_type == "Constant":
            # compute component
            vpsd_comp = c0

        elif component_type == "Harvey":
            # compute component
            vpsd_comp = c0 / (1 + (c1 * freq_avg) ** c2)

        elif component_type == "Lorentz":
            # compute component
            vpsd_comp = c0 * c1**2 / (c1**2 + (freq_avg - c2) ** 2)

        # add component to sum
        vpsd_tot += vpsd_comp

    # logarithmic residuals
    return np.log10(vpsd_avg) - np.log10(vpsd_tot)
