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
    kinematic_bragg_deg,
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
