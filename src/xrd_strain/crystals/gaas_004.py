"""GaAs (004) dynamical diffraction at 10 keV."""

import numpy as np

from xrd_strain.crystals.base import register_crystal
from xrd_strain.crystals.gaas_004_dynamical import xrd_slab_gaas


class GaAs004Calculator:
    name = "gaas_004_10kev"
    substrate_material = "GaAs"

    def compute_intensity(
        self, th_deg: np.ndarray, strain: np.ndarray, dz_angstrom: float, eps: float
    ) -> np.ndarray:
        return xrd_slab_gaas(th_deg, strain, dz_angstrom, eps)


register_crystal(GaAs004Calculator())
