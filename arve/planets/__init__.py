from .detection_test     import detection_test
from .fit_keplerians     import fit_keplerians
from .injection_recovery import injection_recovery
from .plot_detections    import plot_detections
from .plot_keplerians    import plot_keplerians

from typing import Optional

class Planets(
    detection_test,
    fit_keplerians,
    injection_recovery,
    plot_detections,
    plot_keplerians
    ):
    """ARVE Planets sub-class.
    """
    
    def __init__(self, arve):
        self.arve                         = arve
        self.periodograms: Optional[dict] = None
        self.keplerians  : Optional[dict] = None
        self.detections  : Optional[dict] = None