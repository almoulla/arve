from .add_spec         import add_spec
from .add_vrad         import add_vrad
from .compute_vrad_ccf import compute_vrad_ccf

class _Data_classes(add_spec,
                    add_vrad,
                    compute_vrad_ccf):
    pass