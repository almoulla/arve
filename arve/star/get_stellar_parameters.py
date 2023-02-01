"""
get stellar parameters
"""

from   astroquery.simbad import Simbad
import numpy             as     np


def get_stellar_parameters(star):

    # Sun
    if star.target == "Sun":

        star.stellar_parameters["sptype"] = "G2"
        star.stellar_parameters["Teff"  ] = 5770
        star.stellar_parameters["logg"  ] = 4.4
        star.stellar_parameters["Fe_H"  ] = 0.0
        star.stellar_parameters["M"     ] = 1.0
        star.stellar_parameters["R"     ] = 1.0
        star.stellar_parameters["vsini" ] = 1.63

    # other stars
    else:

        # get spectral type
        simbad = Simbad()
        simbad.add_votable_fields("sptype")
        star.stellar_parameters["sptype"] = simbad.query_object(star.target)["SP_TYPE"][0][:2]

        # convert spectral type to number
        sptype_num       =  sptype_to_num(star.stellar_parameters["sptype"])
        sptype_num_table = [sptype_to_num(sptype) for sptype in table["sptype"]]

        # interpolate stellar parameters from table
        keys = ["Teff", "logg", "Fe_H", "M", "R", "vsini"]
        for key in keys:
            star.stellar_parameters[key] = np.interp(sptype_num, sptype_num_table, table[key])
        
    # compute micro- and macro-turbulence
    star.stellar_parameters["vmic"] = 0.85
    star.stellar_parameters["vmac"] = max(0.00, 3.98-(star.stellar_parameters["Teff"]-5770)/650)

# spectral type to number
def sptype_to_num(sptype):
    
    return "OBAFGKM".index(sptype[0])*10 + int(sptype[1])

# table with spectral parameters for main sequence stars
table = np.array(
[
("A0", 9572, 4.3, 0.0, 2.34, 1.80, 255.0),
("A2", 8985, 4.3, 0.0, 2.21, 1.75, 244.0),
("A5", 8306, 4.2, 0.0, 2.04, 1.69, 225.0),
("A7", 7935, 4.2, 0.0, 1.93, 1.68, 210.0),
("F0", 7178, 4.3, 0.0, 1.66, 1.62, 180.0),
("F2", 6909, 4.3, 0.0, 1.56, 1.48, 135.0),
("F5", 6528, 4.3, 0.0, 1.41, 1.40,  20.0),
("F8", 6160, 4.4, 0.0, 1.25, 1.20,   9.0),
("G0", 5943, 4.4, 0.0, 1.16, 1.12,   6.4),
("G2", 5811, 4.4, 0.0, 1.11, 1.08,   4.8),
("G5", 5657, 4.5, 0.0, 1.05, 0.95,   3.4),
("G8", 5486, 4.5, 0.0, 0.97, 0.91,   2.6),
("K0", 5282, 4.6, 0.0, 0.90, 0.83,   2.2),
("K2", 5055, 4.6, 0.0, 0.81, 0.75,   2.0),
("K3", 4973, 4.6, 0.0, 0.79, 0.73,   2.0),
("K5", 4623, 4.6, 0.0, 0.65, 0.64,   1.9),
("K7", 4380, 4.7, 0.0, 0.54, 0.54,   1.7),
("M0", 4212, 4.7, 0.0, 0.46, 0.48,   1.5),
("M2", 4076, 4.7, 0.0, 0.40, 0.43,   0.0),
("M5", 3923, 4.8, 0.0, 0.34, 0.38,   0.0)
],
dtype=[("sptype","U2"), ("Teff","f4"), ("logg","f4"), ("Fe_H","f4"), ("M","f4"), ("R","f4"), ("vsini","f4")]
)