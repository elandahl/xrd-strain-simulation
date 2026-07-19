#!/usr/bin/env python3
"""Cross-code benchmark for Ge/InSb (004) against Stepanov GID_sl.

All curves are raw dynamical reflectivity at 10 keV, sigma polarization,
symmetric Bragg geometry.  No angular or temporal instrument response is
applied.  Cached GID_sl/X0h-database curves and exact input echoes live in
``tests/data``.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from xrd_strain.crystals.base import get_crystal
from xrd_strain.crystals.ge_004_10kev_300k import ge_004_10kev_300k_constants
from xrd_strain.crystals.insb_004_10kev_300k import (
    insb_004_10kev_300k_constants,
)

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "tests" / "data"
REPORT = REPO / "docs" / "gid_sl_ge_insb_benchmark.json"
FIGURE = REPO / "docs" / "images" / "gid_sl_ge_insb_benchmark.png"
DZ_ANGSTROM = 26.7
N_CELLS = 3000
EPS = 1e-6

MATERIALS = {
    "Ge": ("ge", "ge_004_10kev", ge_004_10kev_300k_constants()),
    "InSb": ("insb", "insb_004_10kev", insb_004_10kev_300k_constants()),
}
CASES = {
    "perfect": [],
    "uniform_layer": [(100, 1e-3)],
    "two_step": [(40, 2e-3), (60, 1e-3)],
}


def build_strain(layers: list[tuple[int, float]]) -> np.ndarray:
    strain = np.zeros(N_CELLS)
    start = 0
    for n_cells, value in layers:
        strain[start : start + n_cells] = value
        start += n_cells
    if strain[0] == 0:
        strain[0] = EPS
    return strain


def metrics(scan: np.ndarray, reference: np.ndarray, curve: np.ndarray) -> dict:
    mask = (reference > 1e-6) & (curve > 1e-6)
    dlog = np.log10(curve[mask]) - np.log10(reference[mask])
    return {
        "log10_rms_diff": float(np.sqrt(np.mean(dlog**2))),
        "log10_max_abs_diff": float(np.max(np.abs(dlog))),
        "log10_correlation": float(
            np.corrcoef(np.log10(reference[mask]), np.log10(curve[mask]))[0, 1]
        ),
        "production_peak_arcsec": float(scan[np.argmax(curve)]),
        "gid_sl_peak_arcsec": float(scan[np.argmax(reference)]),
        "production_peak_reflectivity": float(curve.max()),
        "gid_sl_peak_reflectivity": float(reference.max()),
    }


def main() -> int:
    report = {
        "reference": (
            "Stepanov GID_sl/X0h database, (004), 10 keV, sigma, "
            "symmetric Bragg, no instrument"
        ),
        "materials": {},
    }
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), height_ratios=[3, 1])
    all_ok = True

    for col, (material, (file_key, calculator, constants)) in enumerate(
        MATERIALS.items()
    ):
        theta_b = np.degrees(
            np.arcsin(
                constants.wavelength_angstrom
                / (2.0 * (constants.lattice_angstrom / 4.0))
            )
        )
        ax = axes[0, col]
        axr = axes[1, col]
        material_report = {"calculator": calculator, "cases": {}}

        for case, layers in CASES.items():
            scan, gid = np.loadtxt(DATA / f"gid_sl_{file_key}_{case}.dat").T
            strain = build_strain(layers)
            ours = get_crystal(calculator).compute_intensity(
                theta_b + scan / 3600.0, strain, DZ_ANGSTROM, EPS
            )
            result = metrics(scan, gid, ours)
            limit = 0.02 if material == "Ge" else 0.05
            result["passed"] = bool(
                result["log10_rms_diff"] < limit
                and result["log10_correlation"] > 0.999
                and abs(
                    result["production_peak_arcsec"]
                    - result["gid_sl_peak_arcsec"]
                )
                <= 0.5
            )
            all_ok &= result["passed"]
            material_report["cases"][case] = result

            (line,) = ax.semilogy(scan, gid, lw=1.4, label=f"GID_sl: {case}")
            ax.semilogy(
                scan,
                ours,
                ls="--",
                lw=1.0,
                color=line.get_color(),
                label=f"production: {case}",
            )
            mask = gid > 1e-6
            axr.plot(
                scan[mask],
                np.log10(ours[mask]) - np.log10(gid[mask]),
                lw=0.8,
                color=line.get_color(),
                label=case,
            )

        report["materials"][material] = material_report
        ax.set_title(f"{material} (004): solid GID_sl, dashed production")
        ax.set_ylabel("reflectivity")
        ax.set_ylim(1e-7, 2)
        ax.legend(fontsize=7, ncol=2)
        ax.grid(alpha=0.3)
        axr.axhline(0, color="k", lw=0.5)
        axr.set_xlabel(r"$\theta-\theta_B$ (arcsec)")
        axr.set_ylabel(r"$\Delta\log_{10}R$")
        axr.grid(alpha=0.3)

    report["passed"] = all_ok
    fig.suptitle(
        "Ge and InSb production XRD vs Stepanov GID_sl/X0h "
        "(004, 10 keV, sigma, no instrument)"
    )
    fig.tight_layout()
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE, dpi=150)
    plt.close(fig)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")

    for material, mrep in report["materials"].items():
        for case, result in mrep["cases"].items():
            status = "PASS" if result["passed"] else "FAIL"
            print(
                f"[{status}] {material} {case}: "
                f"RMS={result['log10_rms_diff']:.4f}, "
                f"corr={result['log10_correlation']:.6f}, "
                f"peak={result['production_peak_arcsec']:+.2f}\" "
                f"vs {result['gid_sl_peak_arcsec']:+.2f}\""
            )
    print(f"report: {REPORT}")
    print(f"figure: {FIGURE}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
