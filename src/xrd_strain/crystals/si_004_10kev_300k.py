"""Audited Si (004), 10 keV, 300 K structure factors.

Diamond-structure silicon: eight identical basis atoms are in phase for
reflections with h+k+l = 4n, so

* F0 = 8 (Z + f' + i f'')
* Fh = 8 (f0(Q_h) + f' + i f'') DW

Parallel to ``gaas_004_10kev_300k``; GaAs remains the package default.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

ENERGY_KEV = 10.0
WAVELENGTH_ANGSTROM = 1.239841984
# Room-temperature Si lattice constant (Å); standard CODATA-era value.
A_SI_300K = 5.4309
CLASSICAL_ELECTRON_RADIUS_ANGSTROM = 2.8179403205e-5

# Waasmaier & Kirfel (1995) five-Gaussian fit for Si; c chosen so a.sum()+c = Z.
_WK_A = np.array([6.2915, 3.0353, 1.9891, 1.5410, 1.1407])
_WK_B = np.array([2.4386, 32.3337, 0.6785, 81.6937, 1.5496])
_WK_C = 14.0 - float(_WK_A.sum())

# Henke/CXRO at 10,000 eV (interpolated from si.nff): f1≈14.194, f2≈0.214.
FP_SI = 0.1940
FPP_SI = 0.2136

# Sears & Shelley, Acta Cryst. A47, 441 (1991), room-temperature B (Å^2).
B_SI_300K = 0.4613

Z_SI = 14.0


def _f0_wk(s: float) -> float:
    return float(_WK_C + np.sum(_WK_A * np.exp(-_WK_B * s**2)))


@dataclass(frozen=True)
class Si004Constants:
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


def si_004_10kev_300k_constants() -> Si004Constants:
    """Return audited constants for the Si (004) production calculator."""
    a = A_SI_300K
    wavelength = WAVELENGTH_ANGSTROM
    s = 1.0 / (2.0 * (a / 4.0))
    dw = float(np.exp(-B_SI_300K * s**2))
    f0_qh = _f0_wk(s)

    atomic_0 = Z_SI + FP_SI + 1j * FPP_SI
    atomic_h = (f0_qh + FP_SI + 1j * FPP_SI) * dw

    f_0 = 8.0 * atomic_0
    f_h = 8.0 * atomic_h
    f_minus_h = f_h

    scale = -CLASSICAL_ELECTRON_RADIUS_ANGSTROM * wavelength**2 / (np.pi * a**3)

    def susceptibility(factor: complex) -> complex:
        return scale * factor.real - 1j * scale * factor.imag

    return Si004Constants(
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
        f0_qh=f0_qh,
    )
