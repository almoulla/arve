"""
arve structure
"""

from .data.add_rv                import add_rv
from .data.add_spectrum          import add_spectrum

from .star.add_vpsd_component    import add_vpsd_component
from .star.compute_vpsd          import compute_vpsd
from .star.fit_vpsd_coefficients import fit_vpsd_coefficients
from .star.plot_vpsd_components  import plot_vpsd_components

from .planets.add_planet         import add_planet


class ARVE_Structure:
    def __init__(self):
        self.data = _Data(self)
        self.star = _Star(self)
        self.planets = _Planets(self)


class _Data:
    def __init__(self, arve):
        self.arve = arve
        self.rv: dict = {}
        self.spectrum: dict = {}
        self.lightcurve: dict = {}

    def add_rv(self, **kwargs):
        add_rv(self, **kwargs)

    def add_spectrum(self, **kwargs):
        add_spectrum(self, **kwargs)


class _Star:
    def __init__(self, arve):
        self.arve = arve
        self.stellar_parameters: dict = {}
        self.activity_parameters: dict = {}
        self.vpsd: dict = {}
        self.vpsd_components: dict = {}

    def add_vpsd_component(self, **kwargs):
        add_vpsd_component(self, **kwargs)

    def compute_vpsd(self, **kwargs):
        compute_vpsd(self, **kwargs)

    def fit_vpsd_coefficients(self, **kwargs):
        fit_vpsd_coefficients(self, **kwargs)

    def plot_vpsd_components(self, **kwargs):
        plot_vpsd_components(self, **kwargs)


class _Planets:
    def __init__(self, arve):
        self.arve = arve
        self.parameters: dict = {}

    def add_planet(self, **kwargs):
        add_planet(self, **kwargs)