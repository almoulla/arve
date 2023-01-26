"""
compute vpsd
"""


def compute_vpsd(star):

    # RV time series
    time, vrad_val, vrad_err, time_unit, vrad_unit = [star.arve.data.vrad[var] for var in ["time", "vrad_val", "vrad_err", "time_unit", "vrad_unit"]]

    # compute power spectrum
    f, ps, phi, win_f, win_ps, win_area = star.arve.functions.gls_periodogram(time=time, y=vrad_val, sig=vrad_err, normalize=False, win_func=True)

    # compute VPSD
    psd = ps / win_area

    # save in ARVE structure
    star.vpsd = {"freq": f, "vps": ps, "vpsd": psd, "phase": phi}