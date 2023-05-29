from typing import TypeVar

import numpy as np

from arve import ARVE

TPlanets = TypeVar("TPlanets", bound="Planets")


class Planets:
    """ARVE Planets base-class."""

    def __init__(self: TPlanets, arve: ARVE) -> None:
        self.arve = arve
        self.parameters: dict[str, np.float64] = {}

    def add_planet(self: TPlanets) -> None:
        # TODO

        return None
