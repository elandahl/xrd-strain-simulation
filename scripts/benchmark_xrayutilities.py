#!/usr/bin/env python3
"""Tier-4 cross-code benchmark: production calculators vs xrayutilities.

Second *independent* dynamical-diffraction implementation (after Stepanov
GID_sl / X0h, see GID_SL_BENCHMARK.md).  xrayutilities' ``simpack``
``DynamicalModel`` is a generalized 2-beam (4-tiepoint) solver with its own
scattering database, developed and maintained independently of both this
package and the X-ray Server.  Agreement between three independent codes is
much stronger evidence than any pairwise check.

For each material (GaAs 004 and Si 004, 10 keV, sigma polarization, symmetric
Bragg) and each synthetic strain profile we compute three raw rocking curves
(**no** instrument convolution -- ``resolution_width=0`` in xrayutilities and
``instrument=none`` here, because we are comparing diffraction physics, not
the beamline response):

* ``xu``         : xrayutilities DynamicalModel;
* ``production`` : our audited calculator (Waasmaier/Henke/Debye-Waller);
* ``matched``    : our *solver* fed xrayutilities' own chi0/chih.

``production`` vs ``xu`` measures total agreement (implementation + database).
``matched`` vs ``xu`` isolates the numerical dynamical-diffraction engine from
the scattering-database choice -- the same technique used in
FIG3_GID_SL_BENCHMARK.md.

Outputs:
* docs/xrayutilities_benchmark.json    - machine-readable metrics
* docs/images/xrayutilities_benchmark.png - overlays + residuals
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xrayutilities as xu
from xrayutilities.simpack import DynamicalModel, Layer, LayerStack

from xrd_strain.crystals.base import get_crystal
from xrd_strain.crystals.gaas_004_10kev_300k import gaas_004_10kev_300k_constants
from xrd_strain.crystals.gaas_004_dynamical import (
    xrd_slab_gaas_lowmem_with_constants,
)
from xrd_strain.crystals.ge_004_10kev_300k import ge_004_10kev_300k_constants
from xrd_strain.crystals.insb_004_10kev_300k import (
    insb_004_10kev_300k_constants,
)
from xrd_strain.crystals.si_004_10kev_300k import si_004_10kev_300k_constants

REPO = Path(__file__).resolve().parents[1]
REPORT = REPO / "docs" / "xrayutilities_benchmark.json"
FIGURE = REPO / "docs" / "images" / "xrayutilities_benchmark.png"

E_EV = 10000.0
HKL = (0, 0, 4)
DZ_ANGSTROM = 26.7
N_TOTAL_CELLS = 3000
EPS = 1e-6
SCAN = np.linspace(-250.0, 50.0, 601)  # arcsec relative to kinematic Bragg

# (cells, strain) layers, surface first (matches our strain-array convention).
CASES = {
    "perfect": [],
    "uniform_layer": [(100, 1e-3)],
    "two_step": [(40, 2e-3), (60, 1e-3)],
}

MATERIALS = {
    "GaAs": {
        "xu_material": xu.materials.GaAs,
        "constants": gaas_004_10kev_300k_constants(),
        "calculator": "gaas_004_10kev",
    },
    "Si": {
        "xu_material": xu.materials.Si,
        "constants": si_004_10kev_300k_constants(),
        "calculator": "si_004_10kev",
    },
    "Ge": {
        "xu_material": xu.materials.Ge,
        "constants": ge_004_10kev_300k_constants(),
        "calculator": "ge_004_10kev",
    },
    "InSb": {
        "xu_material": xu.materials.InSb,
        "constants": insb_004_10kev_300k_constants(),
        "calculator": "insb_004_10kev",
    },
}


def kinematic_bragg_deg(lattice_a: float, wavelength: float) -> float:
    return float(np.degrees(np.arcsin(wavelength / (2.0 * (lattice_a / 4.0)))))


def factors_from_susceptibilities(
    chi_0: complex, chi_h: complex, constants
) -> np.ndarray:
    """Convert standard chi0/chih into this solver's structure factors."""
    scale = (
        -constants.classical_electron_radius_angstrom
        * constants.wavelength_angstrom**2
        / (np.pi * constants.lattice_angstrom**3)
    )

    def factor(chi: complex) -> complex:
        return chi.real / scale + 1j * (-chi.imag / scale)

    return np.array([factor(chi_0), factor(chi_h)], dtype=np.complex128)


def build_strain(layers: list[tuple[int, float]]) -> np.ndarray:
    strain = np.zeros(N_TOTAL_CELLS)
    start = 0
    for n_cells, value in layers:
        strain[start : start + n_cells] = value
        start += n_cells
    if strain[0] == 0.0:
        strain[0] = EPS  # tiny seed so the perfect crystal is non-degenerate
    return strain


def xu_curve(xu_material, layers: list[tuple[int, float]]) -> tuple[np.ndarray, dict]:
    lam = xu.utilities.en2lam(E_EV)
    theta_b = kinematic_bragg_deg(xu_material.a, lam)
    ai = theta_b + SCAN / 3600.0

    stack = [Layer(xu_material, 1e7)]  # first layer = semi-infinite substrate
    for n_cells, eps in reversed(layers):  # deepest strained layer first
        mat = copy.deepcopy(xu_material)
        strain_matrix = np.zeros((3, 3))
        strain_matrix[2, 2] = eps
        mat.ApplyStrain(strain_matrix)
        stack.append(Layer(mat, n_cells * DZ_ANGSTROM))

    model = DynamicalModel(
        LayerStack("stack", *stack),
        energy=E_EV,
        polarization="S",
        resolution_width=0,
    )
    model.set_hkl(HKL)
    curve = model.simulate(ai, hkl=HKL)
    chi = {
        "chi_0": complex(np.atleast_1d(model.chi0)[0]),
        "chi_h": complex(np.atleast_1d(model.chih["S"])[0]),
    }
    return curve, chi


def our_curve(
    calculator_name: str, constants, factors, layers: list[tuple[int, float]]
) -> np.ndarray:
    theta_b = kinematic_bragg_deg(
        constants.lattice_angstrom, constants.wavelength_angstrom
    )
    th = theta_b + SCAN / 3600.0
    strain = build_strain(layers)
    if factors is None:
        return get_crystal(calculator_name).compute_intensity(
            th, strain, DZ_ANGSTROM, EPS
        )
    return xrd_slab_gaas_lowmem_with_constants(
        th,
        strain,
        DZ_ANGSTROM,
        EPS,
        a_gaas=constants.lattice_angstrom,
        f_gaas=factors,
        re=constants.classical_electron_radius_angstrom,
        wavelength=constants.wavelength_angstrom,
    )


def compare(reference: np.ndarray, curve: np.ndarray) -> dict:
    mask = (reference > 1e-6) & (curve > 1e-6)
    dlog = np.log10(curve[mask]) - np.log10(reference[mask])
    return {
        "log10_rms_diff": float(np.sqrt(np.mean(dlog**2))),
        "log10_max_abs_diff": float(np.max(np.abs(dlog))),
        "log10_correlation": float(
            np.corrcoef(np.log10(reference[mask]), np.log10(curve[mask]))[0, 1]
        ),
        "peak_arcsec": float(SCAN[np.argmax(curve)]),
        "reference_peak_arcsec": float(SCAN[np.argmax(reference)]),
        "peak_reflectivity": float(curve.max()),
        "reference_peak_reflectivity": float(reference.max()),
    }


def main() -> int:
    report: dict = {
        "reference": f"xrayutilities {xu.__version__} DynamicalModel, "
        "S polarization, resolution_width=0",
        "conditions": "10 keV, (004), symmetric Bragg, no instrument convolution",
        "materials": {},
    }

    fig, axes = plt.subplots(
        2, len(MATERIALS), figsize=(20, 8), height_ratios=[3, 1]
    )
    scan_step = float(SCAN[1] - SCAN[0])
    all_ok = True

    for col, (mname, minfo) in enumerate(MATERIALS.items()):
        constants = minfo["constants"]
        material_report = {
            "calculator": minfo["calculator"],
            "susceptibilities": {},
            "cases": {},
        }
        ax = axes[0, col]
        axr = axes[1, col]

        for name, layers in CASES.items():
            xu_int, chi = xu_curve(minfo["xu_material"], layers)
            prod = our_curve(minfo["calculator"], constants, None, layers)
            matched_factors = factors_from_susceptibilities(
                chi["chi_0"], chi["chi_h"], constants
            )
            matched = our_curve(
                minfo["calculator"], constants, matched_factors, layers
            )

            prod_m = compare(xu_int, prod)
            matched_m = compare(xu_int, matched)
            case_ok = bool(
                prod_m["log10_correlation"] > 0.999
                and abs(prod_m["peak_arcsec"] - prod_m["reference_peak_arcsec"])
                <= 2 * scan_step
                and matched_m["log10_rms_diff"] < 0.02
                and matched_m["log10_correlation"] > 0.9995
            )
            all_ok &= case_ok
            material_report["cases"][name] = {
                "layers_cells_strain": layers,
                "production_vs_xu": prod_m,
                "matched_vs_xu": matched_m,
                "passed": case_ok,
            }

            (line,) = ax.semilogy(SCAN, xu_int, lw=1.4, label=f"xu: {name}")
            ax.semilogy(SCAN, prod, ls="--", lw=1.0, color=line.get_color())
            mask = xu_int > 1e-6
            axr.plot(
                SCAN[mask],
                np.log10(prod[mask]) - np.log10(xu_int[mask]),
                lw=0.8,
                color=line.get_color(),
                label=name,
            )

        # record substrate susceptibilities (perfect-crystal chi)
        _, chi0h = xu_curve(minfo["xu_material"], [])
        material_report["susceptibilities"] = {
            "xu_chi_0": [chi0h["chi_0"].real, chi0h["chi_0"].imag],
            "xu_chi_h": [chi0h["chi_h"].real, chi0h["chi_h"].imag],
            "ours_chi_0": [constants.chi_0.real, constants.chi_0.imag],
            "ours_chi_h": [constants.chi_h.real, constants.chi_h.imag],
        }
        report["materials"][mname] = material_report

        ax.set_title(f"{mname} (004): solid = xrayutilities, dashed = production")
        ax.set_ylabel("reflectivity")
        ax.set_ylim(1e-7, 2)
        ax.legend(fontsize=7, ncol=1)
        ax.grid(alpha=0.3)
        axr.axhline(0, color="k", lw=0.5)
        axr.set_ylim(-0.2, 0.2)
        axr.set_xlabel(r"$\theta-\theta_B^{\mathrm{kin}}$ (arcsec)")
        axr.set_ylabel(r"$\Delta\log_{10}R$")
        axr.grid(alpha=0.3)

    report["passed"] = all_ok
    fig.suptitle(
        "Tier-4 cross-code benchmark: production calculators vs xrayutilities "
        "(10 keV, sigma, no instrument)",
        fontsize=12,
    )
    fig.tight_layout()
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE, dpi=150)
    plt.close(fig)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")

    for mname, mrep in report["materials"].items():
        for name, case in mrep["cases"].items():
            status = "PASS" if case["passed"] else "FAIL"
            p = case["production_vs_xu"]
            m = case["matched_vs_xu"]
            print(
                f"[{status}] {mname} {name}: production corr={p['log10_correlation']:.5f} "
                f"RMS={p['log10_rms_diff']:.4f} | matched corr={m['log10_correlation']:.6f} "
                f"RMS={m['log10_rms_diff']:.4f}"
            )
    print(f"report: {REPORT}")
    print(f"figure: {FIGURE}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
