"""Ge (004) dynamical diffraction calculator at 10 keV."""

import numpy as np

from xrd_strain.crystals.base import register_crystal
from xrd_strain.crystals.gaas_004_dynamical import (
    xrd_slab_gaas_lowmem_with_constants,
)
from xrd_strain.crystals.ge_004_10kev_300k import ge_004_10kev_300k_constants

_CONSTANTS = ge_004_10kev_300k_constants()
_THETA_B_DEG = float(
    np.degrees(
        np.arcsin(
            _CONSTANTS.wavelength_angstrom
            / (2.0 * (_CONSTANTS.lattice_angstrom / 4.0))
        )
    )
)


class Ge004Calculator:
    """Ge (004), 10 keV, 300 K production calculator."""

    name = "ge_004_10kev"
    substrate_material = "Ge"
    default_angle_range_deg = (_THETA_B_DEG - 0.05, _THETA_B_DEG + 0.05)

    def compute_intensity(
        self, th_deg: np.ndarray, strain: np.ndarray, dz_angstrom: float, eps: float
    ) -> np.ndarray:
        constants = ge_004_10kev_300k_constants()
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


register_crystal(Ge004Calculator())
