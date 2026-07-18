"""Smoke tests."""

import numpy as np

from xrd_strain.config import XrdConfig
from xrd_strain.io import StrainProfile
from xrd_strain.pipeline import run_xrd


def _dummy_profile() -> StrainProfile:
    n = 100
    z = np.linspace(0, 1e-6, n)
    strain = np.zeros(n)
    strain[50] = 1e-4
    return StrainProfile(
        format_version="1",
        model="test",
        film_material="Cr",
        substrate_material="GaAs",
        z=z,
        displacement=np.zeros(n),
        strain=strain,
        substrate_strain=strain[10:],
        dz=1e-9,
        n_bin_film=10,
        L_film=180e-9,
        L_sub=1800e-9,
        t_max=1e-12,
        metadata={},
    )


def test_xrd_runs_on_dummy_profile():
    result = run_xrd(_dummy_profile(), config=XrdConfig(n_points=10))
    assert len(result.angle_deg) == 10
    assert len(result.intensity) == 10


def test_lowmem_matches_direct_call():
    from xrd_strain.crystals.gaas_004_dynamical import (
        xrd_slab_gaas,
        xrd_slab_gaas_lowmem,
    )

    profile = _dummy_profile()
    th = np.linspace(25.98, 26.08, 25)
    direct = xrd_slab_gaas(th, profile.substrate_strain, 10.0, 1e-6)
    lowmem = xrd_slab_gaas_lowmem(th, profile.substrate_strain, 10.0, 1e-6)
    np.testing.assert_allclose(lowmem, direct, rtol=1e-12)


def test_instrument_options():
    profile = _dummy_profile()
    raw = run_xrd(profile, config=XrdConfig(n_points=30, instrument="none"))
    aps = run_xrd(
        profile,
        config=XrdConfig(
            n_points=30, instrument="aps_7idc", instrument_fwhm_arcsec=20.0
        ),
    )
    # A wide Gaussian must change the curve; both stay finite.
    assert np.all(np.isfinite(raw.intensity))
    assert np.all(np.isfinite(aps.intensity))
    assert not np.allclose(raw.intensity, aps.intensity)
