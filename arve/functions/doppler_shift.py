def doppler_shift(self, wave:list, v:float) -> list:
    """Doppler shift.

    :param wave: wavelengths
    :type wave: list
    :param v: velocity in km/s
    :type v: float
    :return: Doppler-shifted wavelengths
    :rtype: list
    """

    # vacuum speed of light
    c = 2.99792458e5 # [km/s]

    return wave*((1+v/c)/(1-v/c))**(1/2)