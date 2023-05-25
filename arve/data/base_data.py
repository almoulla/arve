from arve import ARVE
from typing import Optional, Union
import numpy as np


class BaseData:
    """ARVE Data base-class."""

    def __init__(self, arve: ARVE) -> None:
        self.arve = arve
        self.spec: dict[str, Union[np.ndarray, str, None]] = {}
        self.vrad: dict[str, Union[np.ndarray, str, None]] = {}

    def add_spec(
        self,
        time: np.ndarray,
        wave: np.ndarray,
        flux_val: np.ndarray,
        flux_err: Optional[np.ndarray] = None,
        time_unit: Optional[str] = None,
        wave_unit: Optional[str] = None,
        flux_unit: Optional[str] = None,
    ) -> None:
        """Add spectral data.

        :param time: time values
        :type time: list
        :param wave: wavelength values
        :type wave: list
        :param flux_val: flux values
        :type flux_val: list
        :param flux_err: flux errors, defaults to None
        :type flux_err: list, optional
        :param time_unit: time unit, defaults to None
        :type time_unit: str, optional
        :param wave_unit: wavelength unit, defaults to None
        :type time_unit: str, optional
        :param flux_unit: flux unit, defaults to None
        :type flux_unit: str, optional
        :return: None
        :rtype: None
        """
        # add dictionary with spectral data
        self.spec = {
            "time": time,
            "wave": wave,
            "flux_val": flux_val,
            "flux_err": flux_err,
            "time_unit": time_unit,
            "wave_unit": wave_unit,
            "flux_unit": flux_unit,
        }

    def add_vrad(
        self,
        time: np.ndarray,
        vrad_val: np.ndarray,
        vrad_err: Optional[np.ndarray] = None,
        time_unit: Optional[str] = None,
        vrad_unit: Optional[str] = None,
    ) -> None:
        """Add radial velocity data.

        :param time: time values
        :param vrad_val: radial velocity values
        :param vrad_err: radial velocity errors, defaults to None
        :param time_unit: time unit, defaults to None
        :param vrad_unit: radial velocity unit, defaults to None
        :return: None
        """
        # add dictionary with radial velocity data
        self.vrad = {
            "time": time,
            "vrad_val": vrad_val,
            "vrad_err": vrad_err,
            "time_unit": time_unit,
            "vrad_unit": vrad_unit,
        }
