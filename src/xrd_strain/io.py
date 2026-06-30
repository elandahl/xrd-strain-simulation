"""Load strain profiles from strain-wave-simulation."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class StrainProfile:
    format_version: str
    model: str
    film_material: str
    substrate_material: str
    z: np.ndarray
    displacement: np.ndarray
    strain: np.ndarray
    substrate_strain: np.ndarray
    dz: float
    n_bin_film: int
    L_film: float
    L_sub: float
    t_max: float
    metadata: dict = field(default_factory=dict)


def load_strain_profile(path: str | Path) -> StrainProfile:
    path = Path(path)
    with np.load(path, allow_pickle=False) as data:
        return StrainProfile(
            format_version=str(data["format_version"]),
            model=str(data["model"]),
            film_material=str(data["film_material"]),
            substrate_material=str(data["substrate_material"]),
            z=data["z"],
            displacement=data["displacement"],
            strain=data["strain"],
            substrate_strain=data["substrate_strain"],
            dz=float(data["dz"]),
            n_bin_film=int(data["n_bin_film"]),
            L_film=float(data["L_film"]),
            L_sub=float(data["L_sub"]),
            t_max=float(data["t_max"]),
            metadata=json.loads(str(data["metadata_json"])),
        )
