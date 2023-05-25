from typing import Optional
import numpy as np


class add_spec:
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
