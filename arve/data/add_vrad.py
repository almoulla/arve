"""
add vrad
"""


def add_vrad(data, time, vrad_val, vrad_err=None, time_unit=None, vrad_unit=None):

    data.vrad = {"time": time, "vrad_val": vrad_val, "vrad_err": vrad_err, "time_unit": time_unit, "vrad_unit": vrad_unit}