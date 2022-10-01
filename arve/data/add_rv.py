"""
add rv
"""


def add_rv(data, label, time, rv, rv_err=None, time_unit=None, rv_unit=None):

    data.rv[label] = {"time": time, "rv": rv, "rv_err": rv_err, "time_unit": time_unit, "rv_unit": rv_unit}