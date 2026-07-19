"""Crystal calculator registry."""

from __future__ import annotations

from typing import Protocol

import numpy as np


class CrystalCalculator(Protocol):
    """Dynamical diffraction backend for a crystal/reflection/energy."""

    name: str
    substrate_material: str
    default_angle_range_deg: tuple[float, float]

    def compute_intensity(
        self, th_deg: np.ndarray, strain: np.ndarray, dz_angstrom: float, eps: float
    ) -> np.ndarray:
        ...


_CALCULATORS: dict[str, CrystalCalculator] = {}


def register_crystal(calc: CrystalCalculator) -> CrystalCalculator:
    _CALCULATORS[calc.name] = calc
    return calc


def get_crystal(name: str) -> CrystalCalculator:
    if name not in _CALCULATORS:
        available = ", ".join(sorted(_CALCULATORS)) or "(none registered)"
        raise KeyError(f"Unknown crystal {name!r}. Available: {available}")
    return _CALCULATORS[name]


def list_crystals() -> list[str]:
    return sorted(_CALCULATORS)
