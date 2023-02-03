def add_vrad(self, time:list, vrad_val:list, vrad_err:list=None, time_unit:str=None, vrad_unit:str=None) -> None:
    """Add radial velocity data.

    :param time: time values
    :type time: list
    :param vrad_val: radial velocity values
    :type vrad_val: list
    :param vrad_err: radial velocity errors, defaults to None
    :type vrad_err: list, optional
    :param time_unit: time unit, defaults to None
    :type time_unit: str, optional
    :param vrad_unit: radial velocity unit, defaults to None
    :type vrad_unit: str, optional
    :return: None
    :rtype: None
    """

    # add dictionary with radial velocity data
    self.vrad = {"time": time, "vrad_val": vrad_val, "vrad_err": vrad_err, "time_unit": time_unit, "vrad_unit": vrad_unit}
    
    return None