"""
AUTHOR : Khaled Al Moulla
DATE   : 2022-01-01

P00    : Collection of plotting standardizations.
"""

#%%
### MODULES

import matplotlib.pylab  as pylab
import matplotlib.pyplot as plt

#%%
### FIGURE PARAMETERS

# Figure size & fontsizes
plot_params = {'figure.figsize'        : (20, 10),
               'figure.titlesize'      :  25,
               'axes.titlesize'        :  25,
               'axes.labelsize'        :  20,
               'xtick.labelsize'       :  20,
               'ytick.labelsize'       :  20,
               'legend.title_fontsize' :  15,
               'legend.fontsize'       :  15
              }
pylab.rcParams.update(plot_params)

#%%
### FUNCTIONS | DESIGNS

# Universal axis parameters
def axis_params(ax=None):
    
    # Define axis
    if ax == None:
        ax = plt.gca()
    
    # X- & Y-axis
    ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    
    # Z-order
    ax.set_zorder(100)

#%%
### FUNCTIONS | PLOTS

# HARPS parameter vs. time
def timeseries(time, param, err, name, ax=None):
    
    # Define axis
    if ax == None:
        ax = plt.gca()
    
    # Plot
    ax.errorbar(time, param, err, fmt=' ', marker='.', color='k', ecolor='r')
    
    # X-axis
    ax.set_xlabel('BJD $-$ 2,400,000 [d]')
    
    # Y-axis
    ax.set_ylabel(name)
    
    # Axis parameters
    axis_params(ax=ax)
    
    # Return axis
    return ax

# Periodogram w/ logarithmic X-axis showing period
def periodogram_logx(period, power, ax=None):
    
    # Define axis
    if ax == None:
        ax = plt.gca()
    
    # Plot
    ax.semilogx(period, power, '-k')
    
    # X-axis
    ax.set_xlabel(r'$P$ [d]')
    
    # Y-axis
    ax.set_ylabel('Power')
    ax.set_ylim(ymin=0)
    
    # Axis parameters
    axis_params(ax=ax)
    
    # Return axis
    return ax

# Periodogram w/ logarithmic axes; X-axis showing frequency
def periodogram_loglog(freq, power, ax=None):
    
    # Define axis
    if ax == None:
        ax = plt.gca()
    
    # Plot
    ax.loglog(freq, power, '-k')
    
    # X-axis
    ax.set_xlabel(r'$f$ $\left[ \mathrm{d}^{-1} \right]$')
    
    # Y-axis
    ax.set_ylabel('Power')
    
    # Axis parameters
    axis_params(ax=ax)
    
    # Return axis
    return ax