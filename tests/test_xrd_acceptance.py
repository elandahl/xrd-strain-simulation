"""Physics acceptance tests for the GaAs (004) dynamical XRD calculator.

These enforce the internal, physics-based checks documented in docs/VALIDATION.md
(perfect-crystal rocking curve, strained-layer Bragg shift, instrument
convolution sanity) plus a bit-for-bit regression against a golden curve
computed from the archival thermo-elastic-gaas calculator.
"""

from pathlib import Path

import numpy as np

from xrd_strain.acceptance import (
    check_instrument_convolution,
    check_perfect_crystal,
    check_strained_layer_shift,
    fwhm_arcsec,
    kinematic_bragg_deg,
)
from xrd_strain.crystals.base import get_crystal
from xrd_strain.crystals.gaas_004_10kev_300k import (
    gaas_004_10kev_300k_constants,
)
from xrd_strain.crystals.gaas_004_dynamical import (
    xrd_slab_gaas,
    xrd_slab_gaas_lowmem,
    xrd_slab_gaas_with_constants,
)

GOLDEN = Path(__file__).parent / "data" / "gaas004_golden.npz"


def test_perfect_crystal_rocking_curve():
    checks = check_perfect_crystal(n_layers=1000, n_angles=400)
    failed = [c.name for c in checks if not c.passed]
    assert not failed, f"failed perfect-crystal checks: {failed}"


def test_strained_layer_bragg_shift():
    checks = check_strained_layer_shift()
    failed = [(c.name, c.detail) for c in checks if not c.passed]
    assert not failed, f"failed strained-layer checks: {failed}"


def test_instrument_convolution_sanity():
    checks = check_instrument_convolution()
    failed = [(c.name, c.detail) for c in checks if not c.passed]
    assert not failed, f"failed instrument checks: {failed}"


def test_frozen_notebook_regression():
    """Current calculator must reproduce the archival golden curve exactly.

    Golden data was generated from thermo-elastic-gaas (tag paper-v1.0)
    xrd_slab_gaas on a fixed strain profile; see docs/VALIDATION.md.
    """
    with np.load(GOLDEN, allow_pickle=False) as data:
        th = data["th"]
        strain = data["strain"]
        dz_ang = float(data["dz_ang"])
        eps = float(data["eps"])
        golden = data["intensity"]

    current = xrd_slab_gaas(th, strain, dz_ang, eps)
    np.testing.assert_allclose(current, golden, rtol=0.0, atol=1e-12)

    # The low-memory path must reproduce it too.
    lowmem = xrd_slab_gaas_lowmem(th, strain, dz_ang, eps)
    np.testing.assert_allclose(lowmem, golden, rtol=0.0, atol=1e-12)

    # The parameterized Tier-3 benchmark path defaults to the exact archival
    # constants and must therefore preserve the same golden curve.
    parameterized = xrd_slab_gaas_with_constants(th, strain, dz_ang, eps)
    np.testing.assert_allclose(parameterized, golden, rtol=0.0, atol=1e-12)


def test_kinematic_bragg_angle_is_sane():
    # GaAs (004) at 10 keV sits near 26.03 deg.
    assert 25.9 < kinematic_bragg_deg() < 26.2


def test_modern_constants_match_x0h_susceptibilities():
    """Externally audited F0 and Fh reproduce Stepanov X0h."""
    c = gaas_004_10kev_300k_constants()

    # X0h, GaAs 004, 10 keV, sigma, X0h/International Tables database.
    x0h_chi0 = -1.8040e-5 + 1j * 3.6087e-7
    x0h_chih = -1.0360e-5 + 1j * 3.3717e-7

    np.testing.assert_allclose(c.chi_0.real, x0h_chi0.real, rtol=0.01)
    np.testing.assert_allclose(c.chi_0.imag, x0h_chi0.imag, rtol=0.06)
    np.testing.assert_allclose(c.chi_h.real, x0h_chih.real, rtol=0.01)
    np.testing.assert_allclose(c.chi_h.imag, x0h_chih.imag, rtol=0.06)


def test_modern_perfect_crystal_matches_x0h_gid_curve():
    """Production calculator matches the independent absorbing X0h/GID curve."""
    c = gaas_004_10kev_300k_constants()
    theta_b = np.degrees(
        np.arcsin(c.wavelength_angstrom / (2.0 * (c.lattice_angstrom / 4.0)))
    )
    th = np.linspace(theta_b - 0.03, theta_b + 0.03, 1000)
    strain = np.zeros(1500)
    strain[0] = 2e-6

    calculator = get_crystal("gaas_004_10kev")
    intensity = calculator.compute_intensity(th, strain, 26.7, 1e-6)

    # Directly measured from the X0h/GID_sl curve queried 2026-07-18.
    assert abs(fwhm_arcsec(th, intensity) - 5.736857) < 0.03
    assert abs(float(intensity.max()) - 0.97282552) < 0.003
