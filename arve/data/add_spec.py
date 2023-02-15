def add_spec(self, time:list, wave:list, flux_val:list, flux_err:list=None, time_unit:str=None, wave_unit:str=None, flux_unit:str=None) -> None:
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
    self.spec = {"time": time, "wave": wave, "flux_val": flux_val, "flux_err": flux_err, "time_unit": time_unit, "wave_unit": wave_unit, "flux_unit": flux_unit}
    
    return None