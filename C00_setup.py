"""
AUTHOR : Khaled Al Moulla
DATE   : 2022-01-01

C00    : Setup.
"""

#%%
### MODULES

import F00_paths as func_paths

#%%
### MODES

target     = 'Sun'
instrument = 'HARPN'
drsversion = 'new'

#%%
### PATHS

path_dir = func_paths.path_local()
path_in  = path_dir + 'Data/Input/'
path_out = path_dir + 'Data/Output/'
path_fig = path_dir + 'Figures/'
path_pap = path_dir + 'Figures/Paper/'

#%%
### JDB LIM

Sjdb_lim = 59000
Njdb_lim = 58315.5