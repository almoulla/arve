"""
AUTHOR : Khaled Al Moulla
DATE   : 2022-01-01

F00    : Directory paths.
"""

#%% MODULES

import os

#%% FUNCTIONS

"""
Path to directory on local computer.

OUTPUT
path : (str) root path
"""
def path_local():
    
    return os.path.expanduser('~/Documents/Work/UNIGE/1. Project/Year 1/Code/Project I/')

"""
Path to directory on university server.

OUTPUT
path : (str) root path
"""
def path_server():
    
    return os.path.expanduser('/hpcstorage/almoulla/Documents/Work/UNIGE/1_Project/Year_1/Code/Project_I/')