"""Acceptance tests for Ge/InSb (004), 10 keV production calculators."""

from pathlib import Path

import numpy as np
import pytest

from xrd_strain.crystals.base import get_crystal, list_crystals
from xrd_strain.crystals.ge_004_10kev_300k import ge_004_10kev_300k_constants
from xrd_strain.crystals.insb_004_10kev_300k import (
    insb_004_10kev_300k_constants,
)

DATA = Path(__file__).parent / "data"

MATERIALS = [
    (
        "Ge",
        "ge",
        "ge_004_10kev",
        ge_004_10kev_300k_constants,
        (25.9, 26.1),
        0.02,
    ),
    (
        "InSb",
        "insb",
        "insb_004_10kev",
        insb_004_10kev_300k_constants,
        (22.4, 22.6),
        0.05,
    ),
]


def _theta_b(constants) -> float:
    return float(
        np.degrees(
            np.arcsin(
                constants.wavelength_angstrom
                / (2.0 * (constants.lattice_angstrom / 4.0))
            )
        )
    )


@pytest.mark.parametrize(
    "material,file_key,calculator,constants_fn,bragg_range,rms_limit", MATERIALS
)
def test_calculator_registered_and_bragg_angle_sane(
    material, file_key, calculator, constants_fn, bragg_range, rms_limit
):
    assert calculator in list_crystals()
    crystal = get_crystal(calculator)
    assert crystal.substrate_material == material
    theta_b = _theta_b(constants_fn())
    assert bragg_range[0] < theta_b < bragg_range[1]
    assert crystal.default_angle_range_deg[0] < theta_b
    assert crystal.default_angle_range_deg[1] > theta_b


@pytest.mark.parametrize(
    "material,file_key,calculator,constants_fn,bragg_range,rms_limit", MATERIALS
)
def test_perfect_crystal_matches_gid_sl(
    material, file_key, calculator, constants_fn, bragg_range, rms_limit
):
    scan, gid = np.loadtxt(DATA / f"gid_sl_{file_key}_perfect.dat").T
    scan, gid = scan[::4], gid[::4]
    constants = constants_fn()
    strain = np.zeros(1500)
    strain[0] = 1e-6
    ours = get_crystal(calculator).compute_intensity(
        _theta_b(constants) + scan / 3600.0, strain, 26.7, 1e-6
    )
    mask = gid > 1e-6
    dlog = np.log10(ours[mask]) - np.log10(gid[mask])
    assert np.sqrt(np.mean(dlog**2)) < rms_limit
    assert np.corrcoef(np.log10(gid[mask]), np.log10(ours[mask]))[0, 1] > 0.999


@pytest.mark.parametrize(
    "material,file_key,calculator,constants_fn,bragg_range,rms_limit", MATERIALS
)
def test_uniform_layer_matches_gid_sl(
    material, file_key, calculator, constants_fn, bragg_range, rms_limit
):
    scan, gid = np.loadtxt(DATA / f"gid_sl_{file_key}_uniform_layer.dat").T
    scan, gid = scan[::4], gid[::4]
    constants = constants_fn()
    strain = np.zeros(3000)
    strain[:100] = 1e-3
    ours = get_crystal(calculator).compute_intensity(
        _theta_b(constants) + scan / 3600.0, strain, 26.7, 1e-6
    )
    mask = gid > 1e-6
    dlog = np.log10(ours[mask]) - np.log10(gid[mask])
    assert np.sqrt(np.mean(dlog**2)) < rms_limit
    assert np.corrcoef(np.log10(gid[mask]), np.log10(ours[mask]))[0, 1] > 0.999
