"""
AUTHOR : Khaled Al Moulla
DATE   : 2021-01-01

F03    : Power spectral density (PSD) componental functions.
"""

#%% MODULES

#%% FUNCTIONS

"""
Lorentz function for periodic component.

INPUT
f   : (array) frequencies
P   : (array) coefficients

OUTPUT
psd : (array) power spectral density at <f>
"""
def psd_P(f, P):
    
    # P[0]: amplitude                [(km/s)^2/d^-1]
    # P[1]: FWHM                     [d^-1]
    # P[2]: central frequency        [d^-1]
    
    psd = P[0]*P[1]**2/(P[1]**2+(f-P[2])**2)
    
    return psd

"""
Harvey function for granulation component.

INPUT
f   : (array) frequencies
G   : (array) coefficients

OUTPUT
psd : (array) power spectral density at <f>
"""
def psd_G(f, G):
    
    # G[0]: amplitude                [(km/s)^2/d^-1]
    # G[1]: characteristic timescale [d]
    # G[2]: power law slope          [unitless]
    
    psd = G[0]/(1+(G[1]*f)**G[2])
    
    return psd