from typing import Callable, Optional, TypeVar, Union

import numpy as np

from arve import ARVE
from data import compute_vrad_ccf

TData = TypeVar("TData", bound="Data")
RT = TypeVar("RT")


def add_method(cls: TData) -> Callable[[Callable[..., RT]], Callable[..., RT]]:
    def decorator(func: Callable[..., RT]) -> Callable[..., RT]:
        setattr(cls, func.__name__, func)
        return func

    return decorator


@add_method(compute_vrad_ccf)
class Data:
    """ARVE Data base-class."""

    def __init__(self: TData, arve: ARVE) -> None:
        self.arve = arve
        self.spec: dict[str, Union[np.ndarray, str, None]] = {}
        self.vrad: dict[str, Union[np.ndarray, str, None]] = {}

    def add_spec(
        self: TData,
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
        :param wave: wavelength values
        :param flux_val: flux values
        :param flux_err: flux errors, defaults to None
        :param time_unit: time unit, defaults to None
        :param wave_unit: wavelength unit, defaults to None
        :param flux_unit: flux unit, defaults to None
        :return: None
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
        self: TData,
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
