"""Type hints for position/variant objects."""
from typing import TypeVar

from .positions import _Fusion, _Position, _SmallVariant, _Variant

Fusion = TypeVar("Fusion", bound=_Fusion)
Position = TypeVar("Position", bound=_Position)
SmallVariant = TypeVar("SmallVariant", bound=_SmallVariant)
Variant = TypeVar("Variant", bound=_Variant)
