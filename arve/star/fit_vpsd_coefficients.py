"""
fit vpsd coefficients
"""

import numpy as np
from lmfit import Parameters, minimize


def fit_vpsd_coefficients(star, label):

    # read VPSD and units
    freq, vpsd = [star.vpsd[label][var] for var in ["freq", "vpsd"]]

    # log-average VPSD
    freq_bin = 10 ** (np.linspace(np.log10(freq[0]), np.log10(freq[-1]), 51))
    freq_avg = (freq_bin[1:] + freq_bin[:-1]) / 2
    vpsd_avg = np.empty(freq_avg.size)
    for i in range(freq_avg.size):
        vpsd_bin = vpsd[(freq > freq_bin[i]) & (freq < freq_bin[i + 1])]
        if len(vpsd_bin) == 0:
            vpsd_avg[i] = np.nan
        else:
            vpsd_avg[i] = np.mean(vpsd_bin)

    # delete NaN values
    i_delete = np.isnan(vpsd_avg)
    freq_avg = np.delete(freq_avg, np.where(i_delete))
    vpsd_avg = np.delete(vpsd_avg, np.where(i_delete))

    # LMFIT parameters
    params = Parameters()

    # loop components
    for comp in star.vpsd_components[label].keys():

        # component dictionary
        comp_dict = star.vpsd_components[label][comp]

        # type, coefficients, and vary
        type = comp_dict["type"]
        coef = comp_dict["coef"]
        vary = comp_dict["vary"]

        # loop coefficients
        for i in range(len(coef)):

            # add parameters
            params.add(comp + "_" + str(i), value=coef[i], vary=vary[i])

    # fit coefficients
    c = minimize(func_res, params, args=(star, label, freq_avg, vpsd_avg))

    # loop components
    for comp in star.vpsd_components[label].keys():

        # component dictionary
        comp_dict = star.vpsd_components[label][comp]

        # coefficients
        coef = comp_dict["coef"]

        # loop coefficients
        for i in range(len(coef)):

            # update coefficients with fitted values
            coef[i] = c.params[comp + "_" + str(i)].value


def func_res(params, star, label, freq_avg, vpsd_avg):

    # empty array for sum of components
    vpsd_tot = np.zeros(len(freq_avg))

    # loop components
    for comp in star.vpsd_components[label].keys():

        comp_dict = star.vpsd_components[label][comp]

        type = comp_dict["type"]

        # component Constant
        if type == "Constant":
            
            # unpack coefficients
            c0 = params[comp + "_0"]

            # compute component
            vpsd_comp = c0

        # component Lorentz
        if type == "Lorentz":

            # unpack coefficients
            c0 = params[comp + "_0"]
            c1 = params[comp + "_1"]
            c2 = params[comp + "_2"]

            # compute component
            vpsd_comp = c0*c1**2/(c1**2+(freq_avg-c2)**2)

        # component Harvey
        if type == "Harvey":

            # unpack coefficients
            c0 = params[comp + "_0"]
            c1 = params[comp + "_1"]
            c2 = params[comp + "_2"]

            # compute component
            vpsd_comp = c0/(1+(c1*freq_avg)**c2)
    
        # add component to sum
        vpsd_tot += vpsd_comp
    
    # logarithmic residual
    logres = np.log10(vpsd_avg) - np.log10(vpsd_tot)

    # return residual
    return logres