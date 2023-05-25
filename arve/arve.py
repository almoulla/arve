import pickle
import random
from pathlib import Path
from typing import TypeVar

from data import Data
from functions import Functions
from planets import Planets
from star import Star

TARVE = TypeVar("TARVE", bound="ARVE")


class ARVE:
    """ARVE main class."""

    def __init__(self: TARVE) -> None:
        self.id: str = str(random.randint(0, 10**8))
        self.data = Data(self)
        self.functions = Functions(self)
        self.planets = Planets(self)
        self.star = Star(self)


def load(arve: str) -> ARVE:
    """Load ARVE object.

    :param arve: ARVE file to load
    :type arve: str
    :return: loaded ARVE object
    :rtype: ARVE
    """
    with Path(arve).open("rb") as file:
        arve = pickle.load(file)
        if not isinstance(arve, ARVE):
            raise TypeError(f"{arve} is not an ARVE object.")
    return arve


def save(arve: ARVE) -> None:
    """Save ARVE object.

    :param arve: ARVE object to save
    :type arve: ARVE
    :return: None
    :rtype: None
    """
    with Path(f"{arve.id}.arve").open("wb") as file:
        pickle.dump(arve, file)


def delete(arve: ARVE) -> None:
    """Delete ARVE object.

    :param arve: ARVE object to delete
    :type arve: ARVE
    :return: None
    :rtype: None
    """
    del arve
