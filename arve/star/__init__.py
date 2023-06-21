from .add_vpsd_components    import add_vpsd_components
from .compute_vpsd           import compute_vpsd
from .fit_vpsd_coefficients  import fit_vpsd_coefficients
from .get_stellar_parameters import get_stellar_parameters
from .plot_vpsd_components   import plot_vpsd_components

from typing import Optional

class Star(
    add_vpsd_components,
    compute_vpsd,
    fit_vpsd_coefficients,
    get_stellar_parameters,
    plot_vpsd_components
    ):
    """ARVE Star sub-class.
    """

    def __init__(self, arve):
        self.arve                               = arve
        self.target            : Optional[str]  = None
        self.stellar_parameters: Optional[dict] = None
        self.vpsd              : Optional[dict] = None
        self.vpsd_components   : Optional[dict] = None