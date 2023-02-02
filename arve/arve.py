"""
arve class
"""

from .data      import add_vrad

from .functions import gls_periodogram

from .planets   import add_planet

from .star      import add_vpsd_components
from .star      import compute_vpsd
from .star      import fit_vpsd_coefficients
from .star      import get_stellar_parameters
from .star      import plot_vpsd_components

import gc
import pickle


class ARVE:

    def __init__(self):
        self.id: str = None
        self.data = _Data(self)
        self.functions = _Functions(self)
        self.planets = _Planets(self)
        self.star = _Star(self)

def load(arve):
    return \
    pickle.load(open(arve, 'rb'))

def save(arve):
    return \
    pickle.dump(arve, open(arve.id+'.arve', 'wb'))

def delete(arve):
    del arve
    gc.collect()


class _Data:

    def __init__(self, arve):
        self.arve = arve
        self.vrad: dict = {}

    def add_vrad(self, **kwargs):
        return \
        add_vrad(self, **kwargs)


class _Functions:

    def __init__(self, arve):
        self.arve = arve
    
    def gls_periodogram(self, **kwargs):
        return \
        gls_periodogram(self, **kwargs)


class _Planets:
    
    def __init__(self, arve):
        self.arve = arve
        self.parameters: dict = {}

    def add_planet(self, **kwargs):
        return \
        add_planet(self, **kwargs)


class _Star:

    def __init__(self, arve):
        self.arve = arve
        self.target: str = None
        self.stellar_parameters: dict = {}
        self.vpsd: dict = {}
        self.vpsd_components: dict = {}

    def add_vpsd_components(self, **kwargs):
        return \
        add_vpsd_components(self, **kwargs)

    def compute_vpsd(self, **kwargs):
        return \
        compute_vpsd(self, **kwargs)

    def fit_vpsd_coefficients(self, **kwargs):
        return \
        fit_vpsd_coefficients(self, **kwargs)

    def get_stellar_parameters(self, **kwargs):
        return \
        get_stellar_parameters(self, **kwargs)

    def plot_vpsd_components(self, **kwargs):
        return \
        plot_vpsd_components(self, **kwargs)