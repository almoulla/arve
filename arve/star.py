from typing import Optional, TypedDict, TypeVar

import matplotlib.pyplot as plt  # type: ignore
import numpy as np
import numpy.typing as npt
from astroquery.simbad import Simbad  # type: ignore
from lmfit import Parameters, minimize  # type: ignore

from arve import ARVE

TStar = TypeVar("TStar", bound="Star")
RT = TypeVar("RT")


class StellarParameters(TypedDict, total=False):
    Teff: float
    logg: float
    M: float
    R: float
    Fe_H: float
    vsini: float
    vmic: float
    vmac: float
    sptype: str


class StellarParametersInterpolated(TypedDict, total=False):
    Teff: npt.NDArray[np.float64]
    logg: npt.NDArray[np.float64]
    M: npt.NDArray[np.float64]
    R: npt.NDArray[np.float64]
    Fe_H: npt.NDArray[np.float64]
    vsini: npt.NDArray[np.float64]
    vmic: npt.NDArray[np.float64]
    vmac: npt.NDArray[np.float64]
    sptype: str


class Vpsd(TypedDict, total=False):
    freq: npt.NDArray[np.float64]
    vps: npt.NDArray[np.float64]
    vpsd: npt.NDArray[np.float64]
    phase: npt.NDArray[np.float64]
    freq_avg: npt.NDArray[np.float64]
    vps_avg: npt.NDArray[np.float64]
    vpsd_avg: npt.NDArray[np.float64]


class VpsdComponent(TypedDict, total=False):
    func_type: str
    coef_val: list[np.float64]
    coef_err: list[int]
    vary: list[bool]


class Star:
    """ARVE Star base-class."""

    def __init__(self: TStar, arve: ARVE) -> None:
        self.arve = arve
        self.target: Optional[str] = None
        self.stellar_parameters: StellarParameters = {}
        self.stellar_parameters_interpolated: StellarParametersInterpolated = {}
        self.vpsd: Vpsd = {}
        self.vpsd_components: dict[str, VpsdComponent] = {}

    def compute_vpsd(self: TStar) -> None:
        """Compute velocity power spectral density (VPSD).

        :return: None
        :rtype: None
        """
        # read RV time series
        time, vrad_val, vrad_err, time_unit, vrad_unit = (
            self.arve.data.vrad[key]
            for key in ["time", "vrad_val", "vrad_err", "time_unit", "vrad_unit"]
        )

        # compute velocity power spectrum
        if not hasattr(self.arve.functions, "gls_periodogram"):
            raise AttributeError("The function 'gls_periodogram' is not available.")
        (
            freq,
            vps,
            phi,
            win_freq,
            win_vps,
            win_area,
        ) = self.arve.functions.gls_periodogram(
            time=time, val=vrad_val, err=vrad_err, normalize=False, win_func=True
        )

        # compute log-average VPS
        freq_bin = 10 ** (np.linspace(np.log10(freq[0]), np.log10(freq[-1]), 51))
        freq_avg = (freq_bin[1:] + freq_bin[:-1]) / 2
        vps_avg = np.empty(freq_avg.size)
        for i in range(freq_avg.size):
            vps_bin = vps[(freq > freq_bin[i]) & (freq < freq_bin[i + 1])]
            vps_avg[i] = np.nan if len(vps_bin) == 0 else np.mean(vps_bin)
        # delete NaN values
        i_delete = np.isnan(vps_avg)
        freq_avg = np.delete(freq_avg, np.where(i_delete))
        vps_avg = np.delete(vps_avg, np.where(i_delete))

        # compute VPSD
        vpsd = vps / win_area
        vpsd_avg = vps_avg / win_area

        # save VPSD
        self.vpsd = {
            "freq": freq,
            "vps": vps,
            "vpsd": vpsd,
            "phase": phi,
            "freq_avg": freq_avg,
            "vps_avg": vps_avg,
            "vpsd_avg": vpsd_avg,
        }

    def add_vpsd_components(
        self: TStar,
        components: list[str] = [
            "Photon_noise",
            "Oscillations",
            "Granulation",
            "Supergranulation",
        ],
    ) -> None:
        """Add velocity power spectral density (VPSD) components.

        :param components: VPSD components, defaults to ["Photon_noise", "Oscillations", "Granulation", "Supergranulation"]
        :return: None
        """
        # read stellar parameters and VPSD
        Teff = self.stellar_parameters["Teff"]
        logg = self.stellar_parameters["logg"]
        M = self.stellar_parameters["M"]
        R = self.stellar_parameters["R"]

        freq = self.vpsd["freq"]
        self.vpsd["vpsd"]
        freq_avg = self.vpsd["freq_avg"]
        vpsd_avg = self.vpsd["vpsd_avg"]

        # photon noise
        if "Photon_noise" in components:
            # coefficients
            c0 = np.min(vpsd_avg)  # amplitude

            # specifications
            name = "Photon_noise"  # name
            func_type = "Constant"  # function type
            coef_val = [c0]  # coefficient values
            coef_err = [0]  # coefficient errors
            vary = [True]  # vary when fitted

            # check noise is non-negative
            if c0 >= 0:
                # save VPSD component
                self.vpsd_components[name] = {
                    "func_type": func_type,
                    "coef_val": coef_val,
                    "coef_err": coef_err,
                    "vary": vary,
                }

        # oscillations
        if "Oscillations" in components:
            # coefficients
            c2 = (
                M / (R**2 * np.sqrt(Teff / 5777)) * 3.05e-3 * 60 * 60 * 24
            )  # central frequency
            c1 = (
                M ** (1 / 2) * R ** (-3 / 2) * 134.9 * 1e-6 * 2 * 60 * 60 * 24
            )  # peak width
            c0 = vpsd_avg[np.argmin(np.abs(freq_avg - c2))]  # amplitude

            # specifications
            name = "Oscillations"  # name
            func_type = "Lorentz"  # function type
            coef_val = [c0, c1, c2]  # coefficient values
            coef_err = [0, 0, 0]  # coefficient errors
            vary = [True, True, True]  # vary when fitted

            # check central frequency is within frequency range
            if (c2 > min(freq)) & (c2 < max(freq)):
                # save VPSD component
                self.vpsd_components[name] = {
                    "func_type": func_type,
                    "coef_val": coef_val,
                    "coef_err": coef_err,
                    "vary": vary,
                }

        # granulation
        if "Granulation" in components:
            # coefficients
            c2 = 2  # exponent
            c1 = (
                1
                / 24
                * (10**logg / 10 ** (4.4)) ** (7 / 9)
                * (Teff / 5777) ** (23 / 9)
            )  # characteristic timescale
            c0 = vpsd_avg[np.argmin(np.abs(freq_avg - 1 / c1))]  # amplitude

            # specifications
            name = "Granulation"  # name
            func_type = "Harvey"  # function type
            coef_val = [c0, c1, c2]  # coefficient values
            coef_err = [0, 0, 0]  # coefficient errors
            vary = [True, True, False]  # vary when fitted

            # check characteristic frequency is within frequency range
            if (1 / c1 > min(freq)) & (1 / c1 < max(freq)):
                # save VPSD component
                self.vpsd_components[name] = {
                    "func_type": func_type,
                    "coef_val": coef_val,
                    "coef_err": coef_err,
                    "vary": vary,
                }

        # supergranulation
        if "Supergranulation" in components:
            # coefficients
            c2 = 2  # exponent
            c1 = (
                10
                / 24
                * (10**logg / 10 ** (4.4)) ** (7 / 9)
                * (Teff / 5777) ** (23 / 9)
            )  # characteristic timescale
            c0 = vpsd_avg[np.argmin(np.abs(freq_avg - 1 / c1))]  # amplitude

            # specifications
            name = "Supergranulation"  # name
            func_type = "Harvey"  # function type
            coef_val = [c0, c1, c2]  # coefficient values
            coef_err = [0, 0, 0]  # coefficient errors
            vary = [True, True, False]  # vary when fitted

            # check characteristic frequency is within frequency range
            if (1 / c1 > min(freq)) & (1 / c1 < max(freq)):
                # save VPSD component
                self.vpsd_components[name] = {
                    "func_type": func_type,
                    "coef_val": coef_val,
                    "coef_err": coef_err,
                    "vary": vary,
                }

    def fit_vpsd_coefficients(self: TStar) -> None:
        """Fit velocity power spectral density (VPSD) coefficients.

        :return: None
        """
        # read VPSD
        self.vpsd["freq"]
        self.vpsd["vpsd"]
        freq_avg = self.vpsd["freq_avg"]
        vpsd_avg = self.vpsd["vpsd_avg"]

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
        params: Parameters,
        self: TStar,
        freq_avg: npt.NDArray[np.float64],
        vpsd_avg: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        # empty array for sum of components
        vpsd_tot = np.zeros(len(freq_avg))

        # loop components
        for comp in self.vpsd_components:
            # component dictionary
            comp_dict = self.vpsd_components[comp]

            # func type
            func_type = comp_dict["func_type"]

            # unpack coefficients
            c0 = params[f"{comp}_0"]
            c1 = params[f"{comp}_1"]
            c2 = params[f"{comp}_2"]

            # type Constant
            if func_type == "Constant":
                # compute component
                vpsd_comp = c0

            elif func_type == "Harvey":
                # compute component
                vpsd_comp = c0 / (1 + (c1 * freq_avg) ** c2)

            elif func_type == "Lorentz":
                # compute component
                vpsd_comp = c0 * c1**2 / (c1**2 + (freq_avg - c2) ** 2)

            # add component to sum
            vpsd_tot += vpsd_comp

        # logarithmic residuals
        return np.log10(vpsd_avg) - np.log10(vpsd_tot)

    def plot_vpsd_components(self: TStar) -> plt.Figure:
        """Plot velocity power spectral density (VPSD) components.

        :return: figure with the VPSD, its log-average and its modeled components
        """
        # figure
        fig = plt.figure()

        # read VPSD and units
        freq = self.vpsd["freq"]
        vpsd = self.vpsd["vpsd"]
        freq_avg = self.vpsd["freq_avg"]
        vpsd_avg = self.vpsd["vpsd_avg"]

        time_unit = self.arve.data.vrad["time_unit"]
        vrad_unit = self.arve.data.vrad["vrad_unit"]

        # plot VPSD and average VPSD
        plt.loglog(freq, vpsd, ls="-", c="k", alpha=0.5)
        plt.loglog(freq_avg, vpsd_avg, ls="None", marker="o", mec="k", mfc="None")

        # empty array for sum of components
        vpsd_tot = np.zeros(len(freq))

        # loop components
        for comp in self.vpsd_components:
            # component dictionary
            comp_dict = self.vpsd_components[comp]

            # func_type and coefficients
            func_type = comp_dict["func_type"]
            coef_val = comp_dict["coef_val"]

            # unpack coefficients
            c0 = coef_val[0]
            c1 = coef_val[1]
            c2 = coef_val[2]

            if func_type == "Constant":
                # compute component
                const_comp = c0

                # plot component
                plt.axhline(const_comp, ls="--", label=comp)

                vpsd_tot += const_comp

            elif func_type == "Harvey":
                # compute component
                harvey_comp = c0 / (1 + (c1 * freq) ** c2)

                # plot component
                plt.loglog(freq, harvey_comp, ls="--", label=comp)

                vpsd_tot += harvey_comp

            elif func_type == "Lorentz":
                # compute component
                lorentz_comp = c0 * c1**2 / (c1**2 + (freq - c2) ** 2)

                # plot component
                plt.loglog(freq, lorentz_comp, ls="--", label=comp)

                vpsd_tot += lorentz_comp

        # plot component sum
        plt.loglog(freq, vpsd_tot, ls="-", c="k", label="Total")

        # plot limits
        plt.xlim(freq[0], freq[-1])

        # plot labels
        plt.xlabel(f"$f$ [{time_unit}$^-1$]")
        plt.ylabel(f"VPSD [({vrad_unit})$^2$ / {time_unit}$^-1$]")

        # plot legend
        leg = plt.legend(loc="lower left")
        leg.set_zorder(101)

        # plot layout
        plt.tight_layout()

        return fig

    def get_stellar_parameters(self: TStar) -> None:
        """Get spectral type from SIMBAD query and stellar parameters from main-sequence table.

        :return: None
        """
        # Sun
        if self.target == "Sun":
            # save stellar parameters
            self.stellar_parameters["sptype"] = "G2"
            self.stellar_parameters["Teff"] = 5770
            self.stellar_parameters["logg"] = 4.4
            self.stellar_parameters["Fe_H"] = 0.0
            self.stellar_parameters["M"] = 1.0
            self.stellar_parameters["R"] = 1.0
            self.stellar_parameters["vsini"] = 1.63

        # other stars
        else:
            # get spectral type from query
            simbad = Simbad()
            simbad.add_votable_fields("sptype")
            self.stellar_parameters["sptype"] = simbad.query_object(self.target)[
                "SP_TYPE"
            ][0][:2]

            # convert spectral types to numbers
            sptype_num = self.arve.functions.sptype_to_num(
                sptype=self.stellar_parameters["sptype"]
            )
            sptype_num_table = [
                self.arve.functions.sptype_to_num(sptype=sptype)
                for sptype in _table["sptype"]
            ]

            # interpolate and save stellar parameters from table
            self.stellar_parameters_interpolated["Teff"] = np.interp(
                sptype_num, sptype_num_table, _table["Teff"]
            )
            self.stellar_parameters_interpolated["logg"] = np.interp(
                sptype_num, sptype_num_table, _table["logg"]
            )
            self.stellar_parameters_interpolated["Fe_H"] = np.interp(
                sptype_num, sptype_num_table, _table["Fe_H"]
            )
            self.stellar_parameters_interpolated["M"] = np.interp(
                sptype_num, sptype_num_table, _table["M"]
            )
            self.stellar_parameters_interpolated["R"] = np.interp(
                sptype_num, sptype_num_table, _table["R"]
            )
            self.stellar_parameters_interpolated["vsini"] = np.interp(
                sptype_num, sptype_num_table, _table["vsini"]
            )
        # compute and save micro- and macro-turbulence
        self.stellar_parameters["vmic"] = 0.85
        self.stellar_parameters["vmac"] = max(
            0.00, 3.98 - (self.stellar_parameters["Teff"] - 5770) / 650
        )


# table with spectral parameters for main sequence stars
_table = np.array(
    [
        ("A0", 9572, 4.3, 0.0, 2.34, 1.80, 255.0),
        ("A2", 8985, 4.3, 0.0, 2.21, 1.75, 244.0),
        ("A5", 8306, 4.2, 0.0, 2.04, 1.69, 225.0),
        ("A7", 7935, 4.2, 0.0, 1.93, 1.68, 210.0),
        ("F0", 7178, 4.3, 0.0, 1.66, 1.62, 180.0),
        ("F2", 6909, 4.3, 0.0, 1.56, 1.48, 135.0),
        ("F5", 6528, 4.3, 0.0, 1.41, 1.40, 20.0),
        ("F8", 6160, 4.4, 0.0, 1.25, 1.20, 9.0),
        ("G0", 5943, 4.4, 0.0, 1.16, 1.12, 6.4),
        ("G2", 5811, 4.4, 0.0, 1.11, 1.08, 4.8),
        ("G5", 5657, 4.5, 0.0, 1.05, 0.95, 3.4),
        ("G8", 5486, 4.5, 0.0, 0.97, 0.91, 2.6),
        ("K0", 5282, 4.6, 0.0, 0.90, 0.83, 2.2),
        ("K2", 5055, 4.6, 0.0, 0.81, 0.75, 2.0),
        ("K3", 4973, 4.6, 0.0, 0.79, 0.73, 2.0),
        ("K5", 4623, 4.6, 0.0, 0.65, 0.64, 1.9),
        ("K7", 4380, 4.7, 0.0, 0.54, 0.54, 1.7),
        ("M0", 4212, 4.7, 0.0, 0.46, 0.48, 1.5),
        ("M2", 4076, 4.7, 0.0, 0.40, 0.43, 0.0),
        ("M5", 3923, 4.8, 0.0, 0.34, 0.38, 0.0),
    ],
    dtype=[
        ("sptype", "U2"),
        ("Teff", "f4"),
        ("logg", "f4"),
        ("Fe_H", "f4"),
        ("M", "f4"),
        ("R", "f4"),
        ("vsini", "f4"),
    ],
)


def _func_res(
    params: Parameters,
    self: Star,
    freq_avg: npt.NDArray[np.float64],
    vpsd_avg: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    # empty array for sum of components
    vpsd_tot = np.zeros(len(freq_avg))

    # loop components
    for comp in self.vpsd_components:
        # component dictionary
        comp_dict = self.vpsd_components[comp]

        # func type
        func_type = comp_dict["func_type"]

        # unpack coefficients
        c0 = params[f"{comp}_0"]
        c1 = params[f"{comp}_1"]
        c2 = params[f"{comp}_2"]

        # type Constant
        if func_type == "Constant":
            # compute component
            vpsd_comp = c0

        elif func_type == "Harvey":
            # compute component
            vpsd_comp = c0 / (1 + (c1 * freq_avg) ** c2)

        elif func_type == "Lorentz":
            # compute component
            vpsd_comp = c0 * c1**2 / (c1**2 + (freq_avg - c2) ** 2)

        # add component to sum
        vpsd_tot += vpsd_comp

    # logarithmic residuals
    return np.log10(vpsd_avg) - np.log10(vpsd_tot)
