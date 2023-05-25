import numpy as np
from typing import Optional


class add_vrad:
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
