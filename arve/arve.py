import gc
import pickle
from typing import Optional

from data import Data

from functions import Functions
from .planets import _Planets_classes
from .star import _Star_classes


class ARVE:
    """ARVE main class."""

    def __init__(self) -> None:
        self.id: Optional[str] = None
        self.data = Data(self)
        self.functions = Functions(self)
        self.planets = _Planets(self)
        self.star = _Star(self)


def load(arve: str) -> ARVE:
    """Load ARVE object.

    :param arve: ARVE file to load
    :type arve: str
    :return: loaded ARVE object
    :rtype: ARVE
    """
    return pickle.load(open(arve, "rb"))


def save(arve: ARVE) -> None:
    """Save ARVE object.

    :param arve: ARVE object to save
    :type arve: ARVE
    :return: None
    :rtype: None
    """
    return pickle.dump(arve, open(arve.id + ".arve", "wb"))


def delete(arve: ARVE) -> None:
    """Delete ARVE object.

    :param arve: ARVE object to delete
    :type arve: ARVE
    :return: None
    :rtype: None
    """
    del arve
    gc.collect()


class _Planets(_Planets_classes):
    """ARVE _Planets sub-class."""

    def __init__(self, arve) -> None:
        self.arve = arve
        self.parameters: dict = {}


class _Star(_Star_classes):
    """ARVE _Star sub-class."""

    def __init__(self, arve) -> None:
        self.arve = arve
        self.target: str = None
        self.stellar_parameters: dict = {}
        self.vpsd: dict = {}
        self.vpsd_components: dict = {}
