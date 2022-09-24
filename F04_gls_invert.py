"""
AUTHOR : Khaled Al Moulla
DATE   : 2021-01-01

F04    : Inverted Generalised Lomb-Scargle (GLS) periodogram.
"""

#%% MODULES

from numpy import arange, zeros, random, \
                  median,                \
                  pi

import compute_syntethic_vrad3 as compute_syntethic_vrad

#%% FUNCTIONS

"""
Reconstructed signal from GLS at given time points.

INPUT
t                : (array) time points on which to reconstruct the signal
ps               : (array) power spectrum
f         (opt.) : (array) frequencies to sum over
phi       (opt.) : (array) phases
c         (opt.) : (array) floating means
ofac      (opt.) : (int)   over-factorization of periodogram

OUTPUT
signal           : (array) reconstructed signal at <t>
"""
def reconstruct(t, ps, f=0, phi=0, c=0, ofac=1):
    
    # Time span and median step
    T  = t[-1] - t[0]
    dt = median(t[1:] - t[:-1])
    
    # Frequencies
    if type(f) == int:
        df = 1/(T*ofac)
        f  = arange(1/T, 1/(2*dt), df)
    
    # Angular frequencies
    omega = 2*pi*f
    
    # Phases
    if type(phi) == int:
        phi = random.uniform(0, 2*pi, len(f))
    
    # Floating means
    if type(c) == int:
        c = zeros(len(f))
    
    # Reconstructed signal
    signal = compute_syntethic_vrad.compute_syntethic_vrad(t, ps, omega, phi, c)
    
    # Return outputs
    return signal