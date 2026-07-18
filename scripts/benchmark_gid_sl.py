#!/usr/bin/env python3
"""Tier-4 cross-code benchmark: production calculator vs Stepanov GID_sl.

Compares the production GaAs (004) @ 10 keV calculator (``gaas_004_10kev``,
instrument = none) against reference rocking curves computed by Sergey
Stepanov's GID_sl dynamical-diffraction server
(https://x-server.gmca.aps.anl.gov/) on two synthetic strain profiles:

1. Uniform strained surface layer: 267 nm of GaAs with da/a = +1e-3 on a
   semi-infinite GaAs substrate.
2. Two-step profile: 106.8 nm at da/a = +2e-3 over 160.2 nm at da/a = +1e-3
   on the substrate.

The GID_sl reference curves are cached in ``tests/data/gid_sl_*.dat`` together
with the exact server input echoes (``.inp``); the query recipe is documented
in docs/GID_SL_BENCHMARK.md so the references can be regenerated.

Outputs:
* docs/gid_sl_benchmark.json  - machine-readable comparison metrics
* docs/images/gid_sl_benchmark.png - overlay + residual figure
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from xrd_strain.crystals.gaas_004_10kev_300k import (
    gaas_004_10kev_300k_constants,
)
from xrd_strain.crystals.gaas_004_dynamical import (
    xrd_slab_gaas_lowmem_with_constants,
)

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "tests" / "data"
REPORT = REPO / "docs" / "gid_sl_benchmark.json"
FIGURE = REPO / "docs" / "images" / "gid_sl_benchmark.png"

# Grid used for both cases: paper-like layer thickness of 26.7 A per cell and
# a slab deep enough (~8 um) to behave as a semi-infinite substrate.
DZ_ANGSTROM = 26.7
N_TOTAL_CELLS = 3000
EPS = 1e-6

# GID_sl profiles were specified in whole angstroms as multiples of DZ.
CASES = {
    "uniform_layer": {
        "dat": DATA / "gid_sl_uniform_layer.dat",
        "description": "267 nm GaAs, da/a = +1e-3, on GaAs substrate",
        "layers": [(100, 1e-3)],  # (cells, strain)
    },
    "two_step": {
        "dat": DATA / "gid_sl_two_step.dat",
        "description": (
            "106.8 nm at da/a = +2e-3 over 160.2 nm at da/a = +1e-3, "
            "on GaAs substrate"
        ),
        "layers": [(40, 2e-3), (60, 1e-3)],
    },
}


def build_strain(layers: list[tuple[int, float]]) -> np.ndarray:
    strain = np.zeros(N_TOTAL_CELLS)
    start = 0
    for n_cells, value in layers:
        strain[start : start + n_cells] = value
        start += n_cells
    return strain


def production_curve(scan_arcsec: np.ndarray, strain: np.ndarray) -> np.ndarray:
    c = gaas_004_10kev_300k_constants()
    theta_b = np.degrees(
        np.arcsin(c.wavelength_angstrom / (2.0 * (c.lattice_angstrom / 4.0)))
    )
    th = theta_b + scan_arcsec / 3600.0
    return xrd_slab_gaas_lowmem_with_constants(
        th,
        strain,
        DZ_ANGSTROM,
        EPS,
        a_gaas=c.lattice_angstrom,
        f_gaas=c.solver_factors,
        re=c.classical_electron_radius_angstrom,
        wavelength=c.wavelength_angstrom,
    )


def compare(scan: np.ndarray, gid: np.ndarray, ours: np.ndarray) -> dict:
    mask = gid > 1e-6
    dlog = np.log10(ours[mask]) - np.log10(gid[mask])
    step = float(scan[1] - scan[0])

    def peak(y: np.ndarray) -> tuple[float, float]:
        i = int(np.argmax(y))
        return float(scan[i]), float(y[i])

    p_gid, r_gid = peak(gid)
    p_ours, r_ours = peak(ours)
    return {
        "scan_step_arcsec": step,
        "substrate_peak_arcsec": {"gid_sl": p_gid, "production": p_ours},
        "substrate_peak_reflectivity": {"gid_sl": r_gid, "production": r_ours},
        "log10_rms_diff": float(np.sqrt(np.mean(dlog**2))),
        "log10_max_abs_diff": float(np.max(np.abs(dlog))),
        "log10_correlation": float(
            np.corrcoef(np.log10(gid[mask]), np.log10(ours[mask]))[0, 1]
        ),
    }


def main() -> int:
    report: dict = {
        "reference": "Stepanov GID_sl (x-server.gmca.aps.anl.gov), queried 2026-07-18",
        "production_calculator": "gaas_004_10kev, instrument=none",
        "grid": {"dz_angstrom": DZ_ANGSTROM, "n_cells": N_TOTAL_CELLS},
        "cases": {},
    }

    fig, axes = plt.subplots(
        2, len(CASES), figsize=(12, 7), sharex="col", height_ratios=[3, 1]
    )

    all_ok = True
    for col, (name, case) in enumerate(CASES.items()):
        scan, gid = np.loadtxt(case["dat"]).T
        strain = build_strain(case["layers"])
        ours = production_curve(scan, strain)
        metrics = compare(scan, gid, ours)
        # Acceptance thresholds mirrored in tests/test_xrd_acceptance.py.
        metrics["passed"] = bool(
            metrics["log10_rms_diff"] < 0.02
            and abs(
                metrics["substrate_peak_arcsec"]["gid_sl"]
                - metrics["substrate_peak_arcsec"]["production"]
            )
            <= metrics["scan_step_arcsec"]
        )
        all_ok &= metrics["passed"]
        report["cases"][name] = {"description": case["description"], **metrics}

        ax = axes[0, col]
        ax.semilogy(scan, gid, "k-", lw=1.6, label="GID_sl (Stepanov server)")
        ax.semilogy(scan, ours, "r--", lw=1.2, label="production gaas_004_10kev")
        ax.set_title(f"{name}: {case['description']}", fontsize=9)
        ax.set_ylabel("reflectivity")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

        axr = axes[1, col]
        mask = gid > 1e-6
        axr.plot(scan[mask], np.log10(ours[mask]) - np.log10(gid[mask]), "b-", lw=0.8)
        axr.axhline(0, color="k", lw=0.5)
        axr.set_ylim(-0.1, 0.1)
        axr.set_xlabel(r"$\theta-\theta_B$ (arcsec)")
        axr.set_ylabel(r"$\Delta\log_{10}R$")
        axr.grid(alpha=0.3)

    fig.suptitle(
        "Tier-4 cross-code benchmark: production calculator vs GID_sl "
        "(GaAs 004, 10 keV, sigma, no instrument)",
        fontsize=11,
    )
    fig.tight_layout()
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE, dpi=150)
    plt.close(fig)

    REPORT.write_text(json.dumps(report, indent=2) + "\n")

    for name, case in report["cases"].items():
        status = "PASS" if case["passed"] else "FAIL"
        print(
            f"[{status}] {name}: log10 RMS {case['log10_rms_diff']:.4f}, "
            f"max {case['log10_max_abs_diff']:.4f}, "
            f"substrate peak {case['substrate_peak_arcsec']['production']:+.2f}\" "
            f"vs GID_sl {case['substrate_peak_arcsec']['gid_sl']:+.2f}\""
        )
    print(f"report: {REPORT}")
    print(f"figure: {FIGURE}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
