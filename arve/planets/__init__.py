from .fit_keplerians     import fit_keplerians
from .injection_recovery import injection_recovery
from .plot_keplerians    import plot_keplerians
from .plot_recoveries    import plot_recoveries
from .recovery_test      import recovery_test

from typing import Optional

class Planets(
    fit_keplerians,
    injection_recovery,
    plot_keplerians,
    plot_recoveries,
    recovery_test
    ):
    """ARVE Planets sub-class.
    """
    
    def __init__(self, arve):
        self.arve                         = arve
        self.periodograms: Optional[dict] = None
        self.keplerians  : Optional[dict] = None
        self.recoveries  : Optional[dict] = None