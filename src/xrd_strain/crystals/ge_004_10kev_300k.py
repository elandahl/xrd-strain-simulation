"""Audited Ge (004), 10 keV, 300 K structure factors.

Diamond-structure germanium has eight identical basis atoms in phase for
(004).  As for the Si production calculator, the forward and diffracting
factors are kept distinct:

* F0 = 8 (Z + f' + i f'')
* Fh = 8 (f0(Q_h) + f' + i f'') exp(-B s^2)

Sources are documented in ``docs/CONSTANTS_PROVENANCE_GE.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

ENERGY_KEV = 10.0
WAVELENGTH_ANGSTROM = 1.239841984
A_GE_300K = 5.65785
CLASSICAL_ELECTRON_RADIUS_ANGSTROM = 2.8179403205e-5

# Waasmaier & Kirfel (1995), evaluated at s = sin(theta)/lambda for Ge (004).
F0H_GE = 20.461447

# Henke/CXRO tables, linearly interpolated at 10,000 eV.
FP_GE = -2.0122025257
FPP_GE = 0.6219585538

# Experimental room-temperature harmonic factor: 0.543(6) A^2 at 294 K.
# J. Harada and S. Sasaki, J. Phys. Soc. Jpn. 38, 866 (1975).
B_GE_300K = 0.543
Z_GE = 32.0


@dataclass(frozen=True)
class Ge004Constants:
    energy_keV: float
    temperature_K: float
    wavelength_angstrom: float
    lattice_angstrom: float
    classical_electron_radius_angstrom: float
    f_0: complex
    f_h: complex
    f_minus_h: complex
    chi_0: complex
    chi_h: complex
    chi_minus_h: complex
    debye_waller: float
    f0_qh: float

    @property
    def solver_factors(self) -> np.ndarray:
        return np.array([self.f_0, self.f_h], dtype=np.complex128)


def ge_004_10kev_300k_constants() -> Ge004Constants:
    """Return audited constants for the Ge (004) production calculator."""
    a = A_GE_300K
    wavelength = WAVELENGTH_ANGSTROM
    s = 1.0 / (2.0 * (a / 4.0))
    dw = float(np.exp(-B_GE_300K * s**2))

    atomic_0 = Z_GE + FP_GE + 1j * FPP_GE
    atomic_h = (F0H_GE + FP_GE + 1j * FPP_GE) * dw
    f_0 = 8.0 * atomic_0
    f_h = 8.0 * atomic_h
    f_minus_h = f_h

    scale = -CLASSICAL_ELECTRON_RADIUS_ANGSTROM * wavelength**2 / (np.pi * a**3)

    def susceptibility(factor: complex) -> complex:
        return scale * factor.real - 1j * scale * factor.imag

    return Ge004Constants(
        energy_keV=ENERGY_KEV,
        temperature_K=300.0,
        wavelength_angstrom=wavelength,
        lattice_angstrom=a,
        classical_electron_radius_angstrom=CLASSICAL_ELECTRON_RADIUS_ANGSTROM,
        f_0=f_0,
        f_h=f_h,
        f_minus_h=f_minus_h,
        chi_0=susceptibility(f_0),
        chi_h=susceptibility(f_h),
        chi_minus_h=susceptibility(f_minus_h),
        debye_waller=dw,
        f0_qh=F0H_GE,
    )
