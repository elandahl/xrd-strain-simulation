"""GaAs (004) dynamical diffraction at 10 keV."""

import numpy as np

from xrd_strain.crystals.base import register_crystal
from xrd_strain.crystals.gaas_004_dynamical import xrd_slab_gaas_lowmem


class GaAs004Calculator:
    name = "gaas_004_10kev"
    substrate_material = "GaAs"

    def compute_intensity(
        self, th_deg: np.ndarray, strain: np.ndarray, dz_angstrom: float, eps: float
    ) -> np.ndarray:
        # Low-memory angle-by-angle evaluation; numerically identical to the
        # original whole-array call (see gaas_004_dynamical.xrd_slab_gaas_lowmem).
        return xrd_slab_gaas_lowmem(th_deg, strain, dz_angstrom, eps)


register_crystal(GaAs004Calculator())
