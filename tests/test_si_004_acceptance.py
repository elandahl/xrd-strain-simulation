"""Acceptance tests for the Si (004) production calculator."""

from pathlib import Path

import numpy as np

from xrd_strain.acceptance import fwhm_arcsec
from xrd_strain.crystals.base import get_crystal, list_crystals
from xrd_strain.crystals.si_004_10kev_300k import si_004_10kev_300k_constants

DATA = Path(__file__).parent / "data"


def test_si_calculator_is_registered_not_default():
    assert "si_004_10kev" in list_crystals()
    assert get_crystal("si_004_10kev").substrate_material == "Si"


def test_si_bragg_angle_is_sane():
    c = si_004_10kev_300k_constants()
    theta_b = np.degrees(
        np.arcsin(c.wavelength_angstrom / (2.0 * (c.lattice_angstrom / 4.0)))
    )
    assert 27.0 < theta_b < 27.4


def test_si_perfect_crystal_matches_gid_sl():
    c = si_004_10kev_300k_constants()
    theta_b = np.degrees(
        np.arcsin(c.wavelength_angstrom / (2.0 * (c.lattice_angstrom / 4.0)))
    )
    scan, gid = np.loadtxt(DATA / "gid_sl_si_perfect.dat").T
    # Use every 2nd point for runtime.
    scan = scan[::2]
    gid = gid[::2]
    th = theta_b + scan / 3600.0
    strain = np.zeros(1500)
    strain[0] = 1e-6
    ours = get_crystal("si_004_10kev").compute_intensity(th, strain, 26.7, 1e-6)

    assert abs(fwhm_arcsec(th, ours) - 2.6699) < 0.15
    assert abs(float(ours.max()) - 0.97596) < 0.01
    assert abs(scan[np.argmax(ours)] - scan[np.argmax(gid)]) <= (scan[1] - scan[0])


def test_si_uniform_layer_matches_gid_sl():
    c = si_004_10kev_300k_constants()
    theta_b = np.degrees(
        np.arcsin(c.wavelength_angstrom / (2.0 * (c.lattice_angstrom / 4.0)))
    )
    scan, gid = np.loadtxt(DATA / "gid_sl_si_uniform_layer.dat").T
    scan = scan[::4]
    gid = gid[::4]
    th = theta_b + scan / 3600.0
    strain = np.zeros(3000)
    strain[:100] = 1e-3
    ours = get_crystal("si_004_10kev").compute_intensity(th, strain, 26.7, 1e-6)
    mask = gid > 1e-6
    dlog = np.log10(ours[mask]) - np.log10(gid[mask])
    assert np.sqrt(np.mean(dlog**2)) < 0.04
    assert np.corrcoef(np.log10(gid[mask]), np.log10(ours[mask]))[0, 1] > 0.999
