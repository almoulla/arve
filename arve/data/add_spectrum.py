"""
add spectrum
"""


def add_spectrum(data, label, wavelength, flux, flux_err=None, wavelength_unit=None, flux_unit=None):

    data.spectrum[label] = {"wavelength": wavelength, "flux": flux, "flux_err": flux_err, "wavelength_unit": wavelength_unit, "flux_err": flux_unit}