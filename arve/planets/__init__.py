from .add_planet import add_planet

from typing import Optional

class Planets(
    add_planet
    ):
    """ARVE Planets sub-class.
    """
    
    def __init__(self, arve):
        self.arve                       = arve
        self.parameters: Optional[dict] = None