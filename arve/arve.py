from   .data      import _Data_classes
from   .functions import _Functions_classes
from   .planets   import _Planets_classes
from   .star      import _Star_classes

import gc
import pickle

class ARVE:
    """ARVE main class.
    """

    def __init__(self):
        self.id: str = None
        self.data = _Data(self)
        self.functions = _Functions(self)
        self.planets = _Planets(self)
        self.star = _Star(self)

def load(arve:str) -> ARVE:
    """Load ARVE object.

    :param arve: ARVE file to load
    :type arve: str
    :return: loaded ARVE object
    :rtype: ARVE
    """

    return pickle.load(open(arve, 'rb'))

def save(arve:ARVE) -> None:
    """Save ARVE object.

    :param arve: ARVE object to save
    :type arve: ARVE
    :return: None
    :rtype: None
    """

    return pickle.dump(arve, open(arve.id+'.arve', 'wb'))

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

class _Data(_Data_classes):
    """ARVE _Data sub-class.
    """

    def __init__(self, arve):
        self.arve = arve
        self.spec: dict = {}
        self.vrad: dict = {}

class _Functions(_Functions_classes):
    """ARVE _Functions sub-class.
    """

    def __init__(self, arve):
        self.arve = arve

class _Planets(_Planets_classes):
    """ARVE _Planets sub-class.
    """
    
    def __init__(self, arve):
        self.arve = arve
        self.parameters: dict = {}

class _Star(_Star_classes):
    """ARVE _Star sub-class.
    """

    def __init__(self, arve):
        self.arve = arve
        self.target: str = None
        self.stellar_parameters: dict = {}
        self.vpsd: dict = {}
        self.vpsd_components: dict = {}