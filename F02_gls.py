"""
AUTHOR : Xavier Dumusque
EDITOR : Khaled Al Moulla
DATE   : 2021-01-01

F02    : The Generalised Lomb-Scargle (GLS) periodogram.
"""

#%% MODULES

import sys

from numpy import arange, ones, \
                  sum, median,  \
                  pi, sin

if (sys.version_info > (3, 0)):
    import periodogram3     as periodogram
    import periodogram_new3 as periodogram_new
else:
    import periodogram
    import periodogram_new

#%% FUNCTIONS

"""
GLS periodogram for data and spectral window function.

INPUT
t                : (array) time points
y                : (array) data at <t>
sig       (opt.) : (array) uncertainties of <y>
ofac      (opt.) : (int)   over-factorization of periodogram
normalize (opt.) : (bool)  normalize power spectrum
win_func  (opt.) : (bool)  compute window function
version   (opt.) : (str)   'old' (Press 1992) or 'new' (Zechmeister & Kürster 2009)

OUTPUT
f                : (array) frequencies
ps               : (array) power spectrum
phi              : (array) phases
win_f            : (array) window function frequencies    (if win_func=True)
win_ps           : (array) window function power spectrum (if win_func=True)
win_area         : (float) window function area           (if win_func=True)
"""
def func_periodogram(t, y, sig=0, ofac=1, normalize=True, win_func=False, version='new'):
    
    # Weights
    if type(sig) == int:
        N = len(t)
        w = ones(N)/N
    else:
        W = sum(1/sig**2)
        w = 1/(W*sig**2)
    
    # Time span and median step
    T  = t[-1] - t[0]
    dt = median(t[1:] - t[:-1])
    
    # Frequency step and span
    df = 1/(T*ofac)
    f  = arange(1/T, 1/(2*dt), df)
    
    # Periodogram
    if version == 'old':
        pg = periodogram    .periodogram(t, y, w, ofac, dt, 1)
    if version == 'new':
        pg = periodogram_new.periodogram(t, y, w, ofac, dt, 1)
    
    # Extract relevant parameters    
    freq, power, phi, a, b, *_ = pg
    
    # Un-normalize
    if normalize == False:
        power = a**2 + b**2
    
    # Spectral window function
    if win_func:
        
        # Central frequency
        fc = f[int(len(f)/2)]
        
        # Sinusoid with unit amplitude
        win_y = sin(2*pi*fc*t)
        
        # Periodogram for window function
        if version == 'old':
             win_pg = periodogram    .periodogram(t, win_y, w, ofac, dt, 1)
        if version == 'new':
             win_pg = periodogram_new.periodogram(t, win_y, w, ofac, dt, 1)
        
        # Extract relevant parameters for window function
        win_freq, win_power, _, win_a, win_b, *_ = win_pg
        
        # Un-normalize window function
        if normalize == False:
            win_power = win_a**2 + win_b**2
        
        # Recenter window function frequencies
        win_freq -= fc
        
        # Window function area
        win_area  = sum(win_power)*df
    
    # Return outputs
    if win_func == True:
        return freq, power, phi, win_freq, win_power, win_area
    if win_func == False:
        return freq, power, phi