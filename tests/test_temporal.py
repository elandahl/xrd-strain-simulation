"""Tests for the temporal instrument response (delay averaging)."""

import numpy as np
import pytest

from xrd_strain.config import XrdConfig
from xrd_strain.io import StrainProfile
from xrd_strain.pipeline import run_xrd
from xrd_strain.temporal import gaussian_delay_weights, run_xrd_delay_averaged


def _profile(peak_idx: int) -> StrainProfile:
    n = 100
    z = np.linspace(0, 1e-6, n)
    strain = np.zeros(n)
    strain[peak_idx] = 1e-4
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


def test_gaussian_delay_weights_normalized_and_centered():
    times, weights = gaussian_delay_weights(0.34e-9, 90e-12, n_samples=9)
    assert len(times) == len(weights) == 9
    np.testing.assert_allclose(weights.sum(), 1.0)
    # symmetric kernel about the nominal delay
    np.testing.assert_allclose(times.mean(), 0.34e-9)
    np.testing.assert_allclose(weights, weights[::-1])
    assert weights.argmax() == 4


def test_zero_fwhm_degenerates_to_single_delay():
    times, weights = gaussian_delay_weights(1e-9, 0.0)
    np.testing.assert_allclose(times, [1e-9])
    np.testing.assert_allclose(weights, [1.0])


def test_average_of_identical_profiles_equals_single_run():
    config = XrdConfig(n_points=30, instrument="none")
    profiles = [_profile(50)] * 5
    weights = np.full(5, 0.2)
    avg = run_xrd_delay_averaged(profiles, weights, config=config)
    single = run_xrd(profiles[0], config=config)
    np.testing.assert_allclose(avg.intensity, single.intensity)


def test_average_is_weighted_linear_combination():
    config = XrdConfig(n_points=30, instrument="none", log10_intensity=False)
    p1, p2 = _profile(40), _profile(60)
    weights = np.array([0.3, 0.7])
    avg = run_xrd_delay_averaged([p1, p2], weights, config=config)
    i1 = run_xrd(p1, config=config).intensity
    i2 = run_xrd(p2, config=config).intensity
    np.testing.assert_allclose(avg.intensity, 0.3 * i1 + 0.7 * i2)


def test_log10_applied_after_averaging_not_before():
    config = XrdConfig(n_points=30, instrument="none", log10_intensity=True)
    p1, p2 = _profile(40), _profile(60)
    weights = np.array([0.5, 0.5])
    avg = run_xrd_delay_averaged([p1, p2], weights, config=config)
    lin = XrdConfig(n_points=30, instrument="none", log10_intensity=False)
    i1 = run_xrd(p1, config=lin).intensity
    i2 = run_xrd(p2, config=lin).intensity
    np.testing.assert_allclose(avg.intensity, np.log10(0.5 * i1 + 0.5 * i2))


def test_mismatched_lengths_raise():
    with pytest.raises(ValueError):
        run_xrd_delay_averaged([_profile(50)], np.array([0.5, 0.5]))
