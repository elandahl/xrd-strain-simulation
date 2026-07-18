"""Audited GaAs (004), 10 keV, 300 K structure factors.

This module deliberately separates the forward and diffracting structure
factors:

* F0 uses f0(Q=0)=Z and controls refraction/absorption (chi_0).
* Fh uses f0(Q_h), anomalous corrections, and species-specific Debye-Waller
  factors and controls diffraction (chi_h).

The original notebook used F0=Fh. That approximation is preserved only in the
legacy calculator.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

ENERGY_KEV = 10.0
WAVELENGTH_ANGSTROM = 1.239841984
A_GAAS_300K = 5.65325
CLASSICAL_ELECTRON_RADIUS_ANGSTROM = 2.8179403205e-5

# Waasmaier & Kirfel (1995), evaluated at s=sin(theta)/lambda for GaAs (004).
F0H_GA = 19.7161
F0H_AS = 21.0693

# Henke/CXRO at 10,000 eV.
FP_GA = -2.9739
FPP_GA = 0.5364
FP_AS = -1.6047
FPP_AS = 0.7101

# Stevenson, Acta Cryst. A50, 621-632 (1994), room-temperature harmonic B.
B_GA_300K = 0.622
B_AS_300K = 0.483

# Atomic numbers: forward non-anomalous scattering f0(Q=0)=Z.
Z_GA = 31.0
Z_AS = 33.0


@dataclass(frozen=True)
class GaAs004Constants:
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
    debye_waller_ga: float
    debye_waller_as: float

    @property
    def solver_factors(self) -> np.ndarray:
        """Return [F0, effective Fh] for the centrosymmetric scalar solver."""
        # For GaAs 004 all basis atoms are in phase and Fh=F-h for the
        # isotropic harmonic model used here.
        return np.array([self.f_0, self.f_h], dtype=np.complex128)


def gaas_004_10kev_300k_constants() -> GaAs004Constants:
    """Return the externally audited constants for the production calculator."""
    a = A_GAAS_300K
    wavelength = WAVELENGTH_ANGSTROM
    s = 1.0 / (2.0 * (a / 4.0))

    dw_ga = float(np.exp(-B_GA_300K * s**2))
    dw_as = float(np.exp(-B_AS_300K * s**2))

    atomic_ga_0 = Z_GA + FP_GA + 1j * FPP_GA
    atomic_as_0 = Z_AS + FP_AS + 1j * FPP_AS
    atomic_ga_h = (F0H_GA + FP_GA + 1j * FPP_GA) * dw_ga
    atomic_as_h = (F0H_AS + FP_AS + 1j * FPP_AS) * dw_as

    # Zincblende (004): all four Ga and four As basis atoms are in phase.
    f_0 = 4.0 * (atomic_ga_0 + atomic_as_0)
    f_h = 4.0 * (atomic_ga_h + atomic_as_h)
    f_minus_h = f_h

    scale = -CLASSICAL_ELECTRON_RADIUS_ANGSTROM * wavelength**2 / (np.pi * a**3)
    # The production solver uses +i absorption internally, while the standard
    # susceptibility convention is chi = -scale_abs*(Re F - i Im F).
    def susceptibility(factor: complex) -> complex:
        return scale * factor.real - 1j * scale * factor.imag

    return GaAs004Constants(
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
        debye_waller_ga=dw_ga,
        debye_waller_as=dw_as,
    )
