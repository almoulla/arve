from .data      import Data
from .functions import Functions
from .planets   import Planets
from .star      import Star

import gc
import pickle

from typing import Optional

class ARVE:
    """ARVE main class.
    """

    def __init__(self):
        self.id       : Optional[str] = None
        self.data     : Data          = Data(self)
        self.functions: Functions     = Functions(self)
        self.planets  : Planets       = Planets(self)
        self.star     : Star          = Star(self)

def save(arve:ARVE) -> None:
    """Save ARVE object.

    :param arve: ARVE object to save
    :type arve: ARVE
    :return: None
    :rtype: None
    """

    return pickle.dump(arve, open(arve.id+'.arve', 'wb'))

def load(arve:str) -> ARVE:
    """Load ARVE object.

    :param arve: ARVE file to load
    :type arve: str
    :return: loaded ARVE object
    :rtype: ARVE
    """

    return pickle.load(open(arve, 'rb'))

def delete(arve:ARVE) -> None:
    """Delete ARVE object.

    :param arve: ARVE object to delete
    :type arve: ARVE
    :return: None
    :rtype: None
    """
    
    del arve
    gc.collect()
    
    return None