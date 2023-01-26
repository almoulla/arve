"""
arve class
"""

from .data.add_vrad              import add_vrad

from .functions.gls_periodogram  import gls_periodogram

from .planets.add_planet         import add_planet

from .star.add_vpsd_component    import add_vpsd_component
from .star.compute_vpsd          import compute_vpsd
from .star.fit_vpsd_coefficients import fit_vpsd_coefficients
from .star.plot_vpsd_components  import plot_vpsd_components


class ARVE:

    def __init__(self):
        self.data = _Data(self)
        self.functions = _Functions(self)
        self.planets = _Planets(self)
        self.star = _Star(self)


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
        self.stellar_parameters: dict = {}
        self.vpsd: dict = {}
        self.vpsd_components: dict = {}

    def add_vpsd_component(self, **kwargs):
        return \
        add_vpsd_component(self, **kwargs)

    def compute_vpsd(self, **kwargs):
        return \
        compute_vpsd(self, **kwargs)

    def fit_vpsd_coefficients(self, **kwargs):
        return \
        fit_vpsd_coefficients(self, **kwargs)

    def plot_vpsd_components(self, **kwargs):
        return \
        plot_vpsd_components(self, **kwargs)