from .add_data          import add_data
from .compute_spec_mast import compute_spec_mast
from .compute_vrad_ccf  import compute_vrad_ccf
from .compute_vrad_lbl  import compute_vrad_lbl
from .get_aux_data      import get_aux_data
from .plot_vrad         import plot_vrad
from .read_spec         import read_spec

from typing import Optional

class Data(
    add_data,
    compute_spec_mast,
    compute_vrad_ccf,
    compute_vrad_lbl,
    get_aux_data,
    plot_vrad,
    read_spec
    ):
    """ARVE Data sub-class.
    """

    def __init__(self, arve):
        self.arve                     = arve
        self.time    : Optional[dict] = None
        self.vrad    : Optional[dict] = None
        self.spec    : Optional[dict] = None
        self.aux_data: Optional[dict] = None