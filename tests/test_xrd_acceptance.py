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
    xrd_slab_gaas_lowmem_with_constants,
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


def test_production_matches_gid_sl_on_strained_layers():
    """Tier-4 cross-code check: production calculator vs Stepanov GID_sl.

    Reference curves were computed by the GID_sl server (queried 2026-07-18)
    for two synthetic strain profiles on a GaAs substrate; the exact server
    inputs are cached next to the .dat files. See docs/GID_SL_BENCHMARK.md.
    """
    c = gaas_004_10kev_300k_constants()
    theta_b = np.degrees(
        np.arcsin(c.wavelength_angstrom / (2.0 * (c.lattice_angstrom / 4.0)))
    )
    dz = 26.7

    cases = {
        "gid_sl_uniform_layer.dat": [(100, 1e-3)],
        "gid_sl_two_step.dat": [(40, 2e-3), (60, 1e-3)],
    }
    for fname, layers in cases.items():
        scan, gid = np.loadtxt(Path(__file__).parent / "data" / fname).T
        strain = np.zeros(3000)
        start = 0
        for n_cells, value in layers:
            strain[start : start + n_cells] = value
            start += n_cells

        th = theta_b + scan / 3600.0
        ours = xrd_slab_gaas_lowmem_with_constants(
            th,
            strain,
            dz,
            1e-6,
            a_gaas=c.lattice_angstrom,
            f_gaas=c.solver_factors,
            re=c.classical_electron_radius_angstrom,
            wavelength=c.wavelength_angstrom,
        )

        mask = gid > 1e-6
        dlog = np.log10(ours[mask]) - np.log10(gid[mask])
        assert np.sqrt(np.mean(dlog**2)) < 0.02, fname
        # Substrate peak must land on the same scan sample.
        step = scan[1] - scan[0]
        assert abs(scan[np.argmax(ours)] - scan[np.argmax(gid)]) <= step, fname


def test_fig3_strain_matches_gid_sl_when_constants_are_controlled():
    """Realistic d'Alembert strain agrees cross-code with like-for-like chi.

    The production curve is also checked, but its expected residual includes
    the deliberate Waasmaier/Henke/Stevenson vs X0h database difference.
    Every tenth cached scan point is used to keep the test runtime practical.
    """
    data_dir = Path(__file__).parent / "data"
    with np.load(
        data_dir / "fig3_dalembert_substrate_strain.npz", allow_pickle=False
    ) as data:
        strain = data["substrate_strain"]
        dz = float(data["dz_angstrom"])
    scan, gid = np.loadtxt(data_dir / "gid_sl_fig3_dalembert.dat").T
    scan = scan[::10]
    gid = gid[::10]

    c = gaas_004_10kev_300k_constants()
    theta_b = np.degrees(
        np.arcsin(c.wavelength_angstrom / (2.0 * (c.lattice_angstrom / 4.0)))
    )
    theta = theta_b + scan / 3600.0

    scale = (
        -c.classical_electron_radius_angstrom
        * c.wavelength_angstrom**2
        / (np.pi * c.lattice_angstrom**3)
    )
    x0h_chi0 = -1.8040e-5 + 1j * 3.6087e-7
    x0h_chih = -1.0360e-5 + 1j * 3.3717e-7

    def factor(chi):
        return chi.real / scale + 1j * (-chi.imag / scale)

    curves = {}
    for name, factors in {
        "production": c.solver_factors,
        "gid_constants_matched": np.array(
            [factor(x0h_chi0), factor(x0h_chih)], dtype=np.complex128
        ),
    }.items():
        curves[name] = xrd_slab_gaas_lowmem_with_constants(
            theta,
            strain,
            dz,
            1e-6,
            a_gaas=c.lattice_angstrom,
            f_gaas=factors,
            re=c.classical_electron_radius_angstrom,
            wavelength=c.wavelength_angstrom,
        )

    mask = gid > 1e-6
    log_gid = np.log10(gid[mask])
    production_log = np.log10(curves["production"][mask])
    matched_log = np.log10(curves["gid_constants_matched"][mask])
    production_rms = np.sqrt(np.mean((production_log - log_gid) ** 2))
    matched_rms = np.sqrt(np.mean((matched_log - log_gid) ** 2))

    assert np.corrcoef(log_gid, production_log)[0, 1] > 0.994
    assert matched_rms < 0.03
    assert np.corrcoef(log_gid, matched_log)[0, 1] > 0.999
    assert matched_rms < production_rms / 3.0


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
