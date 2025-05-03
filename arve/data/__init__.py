from .add_data               import add_data
from .compute_spec_reference import compute_spec_reference
from .compute_vrad_ccf       import compute_vrad_ccf
from .compute_vrad_lbl       import compute_vrad_lbl
from .get_aux_data           import get_aux_data
from .plot_spec_data         import plot_spec_data
from .plot_vrad              import plot_vrad
from .read_spec              import read_spec

from typing import Optional, Any

class Data(
    add_data,
    compute_spec_reference,
    compute_vrad_ccf,
    compute_vrad_lbl,
    get_aux_data,
    plot_spec_data,
    plot_vrad,
    read_spec
    ):
    """ARVE Data subclass.
    """

    def __init__(self, arve):
        self.arve                                      = arve
        self.aux_data        : Optional[dict[str,Any]] = None
        self.time            : Optional[dict[str,Any]] = None
        self.vrad            : Optional[dict[str,Any]] = None
        self.vrad_components : Optional[dict[str,Any]] = None
        self.spec            : Optional[dict[str,Any]] = None
        self.spec_reference  : Optional[dict[str,Any]] = None
        self.ccf             : Optional[dict[str,Any]] = None