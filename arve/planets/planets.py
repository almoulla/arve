from typing import TypeVar

from arve import ARVE

TPlanets = TypeVar("TPlanets", bound="Planets")


class Planets:
    """ARVE Planets base-class."""

    def __init__(self: TPlanets, arve: ARVE) -> None:
        self.arve = arve
        self.parameters: dict = {}

    def add_planet(self: TPlanets) -> None:
        # TODO

        return None
