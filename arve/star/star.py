from typing import Callable, Optional, Type, TypeVar

from add_vpsd_components import add_vpsd_components
from compute_vpsd import compute_vpsd
from fit_vpsd_coefficients import fit_vpsd_coefficients
from get_stellar_parameters import get_stellar_parameters
from plot_vpsd_components import plot_vpsd_components

from arve import ARVE

TStar = TypeVar("TStar", bound="Star")
RT = TypeVar("RT")


def add_methods(functions: list[Callable[..., RT]]) -> Callable[..., Type["Star"]]:
    """Add methods to the base class."""

    def decorator(cls: Type[Star]) -> Type[Star]:
        for function in functions:
            setattr(cls, function.__name__, function)
        return cls

    return decorator


@add_methods(
    [
        compute_vpsd,
        add_vpsd_components,
        fit_vpsd_coefficients,
        plot_vpsd_components,
        get_stellar_parameters,
    ]
)
class Star:
    """ARVE Star base-class."""

    def __init__(self: TStar, arve: ARVE) -> None:
        self.arve = arve
        self.target: Optional[str] = None
        self.stellar_parameters: dict = {}
        self.vpsd: dict = {}
        self.vpsd_components: dict = {}
