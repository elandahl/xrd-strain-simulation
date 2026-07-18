"""GaAs (004) dynamical diffraction calculators at 10 keV."""

import numpy as np

from xrd_strain.crystals.base import register_crystal
from xrd_strain.crystals.gaas_004_10kev_300k import (
    gaas_004_10kev_300k_constants,
)
from xrd_strain.crystals.gaas_004_dynamical import (
    xrd_slab_gaas_lowmem,
    xrd_slab_gaas_lowmem_with_constants,
)


class GaAs004Calculator:
    """Externally benchmarked GaAs (004), 10 keV, 300 K calculator."""

    name = "gaas_004_10kev"
    substrate_material = "GaAs"

    def compute_intensity(
        self, th_deg: np.ndarray, strain: np.ndarray, dz_angstrom: float, eps: float
    ) -> np.ndarray:
        constants = gaas_004_10kev_300k_constants()
        return xrd_slab_gaas_lowmem_with_constants(
            th_deg,
            strain,
            dz_angstrom,
            eps,
            a_gaas=constants.lattice_angstrom,
            f_gaas=constants.solver_factors,
            re=constants.classical_electron_radius_angstrom,
            wavelength=constants.wavelength_angstrom,
        )


class GaAs004LegacyCalculator:
    """Notebook-faithful calculator retained only for archival reproduction."""

    name = "gaas_004_10kev_legacy"
    substrate_material = "GaAs"

    def compute_intensity(
        self, th_deg: np.ndarray, strain: np.ndarray, dz_angstrom: float, eps: float
    ) -> np.ndarray:
        return xrd_slab_gaas_lowmem(th_deg, strain, dz_angstrom, eps)


register_crystal(GaAs004Calculator())
register_crystal(GaAs004LegacyCalculator())
