"""Cross-code acceptance against xrayutilities (second external code).

Skipped automatically when xrayutilities is not installed (it needs an
OpenMP-capable build; see docs/XU_BENCHMARK.md). When present, it enforces:

* our *numerical* engine reproduces xrayutilities' DynamicalModel to
  <0.01 log10 RMS when fed the same susceptibilities (isolates numerics from
  the scattering database);
* the audited *production* calculator keeps >0.999 log10 correlation with
  xrayutilities across the whole curve (shape), the residual being the
  deliberate Debye-Waller database difference.

Kept small (coarse scan, GaAs + Si perfect and one strained layer) so it runs
in a few seconds.
"""

import copy

import numpy as np
import pytest

xu = pytest.importorskip("xrayutilities")

from xrayutilities.simpack import DynamicalModel, Layer, LayerStack  # noqa: E402

from xrd_strain.crystals.base import get_crystal  # noqa: E402
from xrd_strain.crystals.gaas_004_10kev_300k import (  # noqa: E402
    gaas_004_10kev_300k_constants,
)
from xrd_strain.crystals.gaas_004_dynamical import (  # noqa: E402
    xrd_slab_gaas_lowmem_with_constants,
)
from xrd_strain.crystals.ge_004_10kev_300k import (  # noqa: E402
    ge_004_10kev_300k_constants,
)
from xrd_strain.crystals.insb_004_10kev_300k import (  # noqa: E402
    insb_004_10kev_300k_constants,
)
from xrd_strain.crystals.si_004_10kev_300k import (  # noqa: E402
    si_004_10kev_300k_constants,
)

E_EV = 10000.0
HKL = (0, 0, 4)
DZ = 26.7
SCAN = np.linspace(-200.0, 40.0, 161)  # arcsec, coarse for speed


def _bragg(a, lam):
    return float(np.degrees(np.arcsin(lam / (2.0 * (a / 4.0)))))


def _xu_curve(material, layers):
    lam = xu.utilities.en2lam(E_EV)
    ai = _bragg(material.a, lam) + SCAN / 3600.0
    stack = [Layer(material, 1e7)]
    for n_cells, eps in reversed(layers):
        mat = copy.deepcopy(material)
        S = np.zeros((3, 3))
        S[2, 2] = eps
        mat.ApplyStrain(S)
        stack.append(Layer(mat, n_cells * DZ))
    model = DynamicalModel(
        LayerStack("s", *stack), energy=E_EV, polarization="S", resolution_width=0
    )
    model.set_hkl(HKL)
    curve = model.simulate(ai, hkl=HKL)
    chi0 = complex(np.atleast_1d(model.chi0)[0])
    chih = complex(np.atleast_1d(model.chih["S"])[0])
    return curve, chi0, chih


def _strain(layers, n=2000):
    s = np.zeros(n)
    start = 0
    for n_cells, eps in layers:
        s[start : start + n_cells] = eps
        start += n_cells
    if s[0] == 0:
        s[0] = 1e-6
    return s


def _factors(chi0, chih, c):
    scale = (
        -c.classical_electron_radius_angstrom
        * c.wavelength_angstrom**2
        / (np.pi * c.lattice_angstrom**3)
    )
    return np.array(
        [
            chi0.real / scale - 1j * chi0.imag / scale,
            chih.real / scale - 1j * chih.imag / scale,
        ],
        dtype=np.complex128,
    )


def _corr_rms(ref, curve):
    m = (ref > 1e-6) & (curve > 1e-6)
    dlog = np.log10(curve[m]) - np.log10(ref[m])
    corr = np.corrcoef(np.log10(ref[m]), np.log10(curve[m]))[0, 1]
    return float(corr), float(np.sqrt(np.mean(dlog**2)))


CASES = [
    ("GaAs", xu.materials.GaAs, gaas_004_10kev_300k_constants(), "gaas_004_10kev", []),
    (
        "GaAs",
        xu.materials.GaAs,
        gaas_004_10kev_300k_constants(),
        "gaas_004_10kev",
        [(100, 1e-3)],
    ),
    ("Si", xu.materials.Si, si_004_10kev_300k_constants(), "si_004_10kev", []),
    ("Ge", xu.materials.Ge, ge_004_10kev_300k_constants(), "ge_004_10kev", []),
    (
        "Ge",
        xu.materials.Ge,
        ge_004_10kev_300k_constants(),
        "ge_004_10kev",
        [(100, 1e-3)],
    ),
    (
        "InSb",
        xu.materials.InSb,
        insb_004_10kev_300k_constants(),
        "insb_004_10kev",
        [],
    ),
    (
        "InSb",
        xu.materials.InSb,
        insb_004_10kev_300k_constants(),
        "insb_004_10kev",
        [(100, 1e-3)],
    ),
]


@pytest.mark.parametrize("mname,material,constants,calc,layers", CASES)
def test_matched_susceptibilities_reproduce_xrayutilities(
    mname, material, constants, calc, layers
):
    xu_int, chi0, chih = _xu_curve(material, layers)
    th = _bragg(constants.lattice_angstrom, constants.wavelength_angstrom) + SCAN / 3600.0
    strain = _strain(layers)
    matched = xrd_slab_gaas_lowmem_with_constants(
        th,
        strain,
        DZ,
        1e-6,
        a_gaas=constants.lattice_angstrom,
        f_gaas=_factors(chi0, chih, constants),
        re=constants.classical_electron_radius_angstrom,
        wavelength=constants.wavelength_angstrom,
    )
    corr, rms = _corr_rms(xu_int, matched)
    assert corr > 0.9995, f"{mname}: matched corr {corr}"
    assert rms < 0.01, f"{mname}: matched RMS {rms}"


@pytest.mark.parametrize("mname,material,constants,calc,layers", CASES)
def test_production_shape_matches_xrayutilities(
    mname, material, constants, calc, layers
):
    xu_int, _, _ = _xu_curve(material, layers)
    th = _bragg(constants.lattice_angstrom, constants.wavelength_angstrom) + SCAN / 3600.0
    strain = _strain(layers)
    prod = get_crystal(calc).compute_intensity(th, strain, DZ, 1e-6)
    corr, _ = _corr_rms(xu_int, prod)
    assert corr > 0.999, f"{mname}: production corr {corr}"
