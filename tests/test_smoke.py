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
