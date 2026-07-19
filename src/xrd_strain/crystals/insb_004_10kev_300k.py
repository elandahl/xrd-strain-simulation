"""Audited InSb (004), 10 keV, 300 K structure factors.

InSb is zincblende.  All four In and four Sb basis atoms are in phase for
(004), but the forward and diffracting factors remain distinct:

* F0 = 4 [(Z_In + f'_In + i f''_In) + (Z_Sb + f'_Sb + i f''_Sb)]
* Fh = 4 [(f0_In(Q_h) + f'_In + i f''_In) DW_In
          + (f0_Sb(Q_h) + f'_Sb + i f''_Sb) DW_Sb]

The room-temperature displacement parameters for InSb are less securely
known than for Ge, Si, or GaAs.  The production value below is the effective
isotropic Debye-model B at 300 K for theta_D=160 K; its provenance and
uncertainty are explicit in ``docs/CONSTANTS_PROVENANCE_INSB.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

ENERGY_KEV = 10.0
WAVELENGTH_ANGSTROM = 1.239841984
A_INSB_300K = 6.47937
CLASSICAL_ELECTRON_RADIUS_ANGSTROM = 2.8179403205e-5

# Waasmaier & Kirfel (1995), evaluated at s = sin(theta)/lambda for (004).
F0H_IN = 33.601616
F0H_SB = 35.023630

# Henke/CXRO tables, linearly interpolated at 10,000 eV.
FP_IN = 0.0815480823
FPP_IN = 3.4190647240
FP_SB = 0.0753767072
FPP_SB = 3.9183544060

# Effective isotropic 300 K value from a Debye model with theta_D=160 K and
# the mean In/Sb atomic mass.  Reid & Pirie (Acta Cryst. A39, 1-13, 1983)
# document the material/model dependence of species-resolved InSb factors.
B_IN_300K = 1.14730
B_SB_300K = 1.14730

Z_IN = 49.0
Z_SB = 51.0


@dataclass(frozen=True)
class InSb004Constants:
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
    debye_waller_in: float
    debye_waller_sb: float
    f0_qh_in: float
    f0_qh_sb: float

    @property
    def solver_factors(self) -> np.ndarray:
        return np.array([self.f_0, self.f_h], dtype=np.complex128)


def insb_004_10kev_300k_constants() -> InSb004Constants:
    """Return audited constants for the InSb (004) production calculator."""
    a = A_INSB_300K
    wavelength = WAVELENGTH_ANGSTROM
    s = 1.0 / (2.0 * (a / 4.0))
    dw_in = float(np.exp(-B_IN_300K * s**2))
    dw_sb = float(np.exp(-B_SB_300K * s**2))

    atomic_in_0 = Z_IN + FP_IN + 1j * FPP_IN
    atomic_sb_0 = Z_SB + FP_SB + 1j * FPP_SB
    atomic_in_h = (F0H_IN + FP_IN + 1j * FPP_IN) * dw_in
    atomic_sb_h = (F0H_SB + FP_SB + 1j * FPP_SB) * dw_sb

    f_0 = 4.0 * (atomic_in_0 + atomic_sb_0)
    f_h = 4.0 * (atomic_in_h + atomic_sb_h)
    f_minus_h = f_h

    scale = -CLASSICAL_ELECTRON_RADIUS_ANGSTROM * wavelength**2 / (np.pi * a**3)

    def susceptibility(factor: complex) -> complex:
        return scale * factor.real - 1j * scale * factor.imag

    return InSb004Constants(
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
        debye_waller_in=dw_in,
        debye_waller_sb=dw_sb,
        f0_qh_in=F0H_IN,
        f0_qh_sb=F0H_SB,
    )
