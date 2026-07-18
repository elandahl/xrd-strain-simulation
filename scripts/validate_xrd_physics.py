#!/usr/bin/env python3
"""Run the internal physics acceptance suite for the GaAs (004) XRD calculator.

These are model-internal, physics-based sanity checks that should be green
*before* any external benchmarking (Sci. Rep. Fig. 3, APS/PLS data):

  * unstrained perfect-crystal rocking curve (Bragg position, reflectivity,
    Darwin width),
  * uniformly strained surface layer -> Bragg peak shift Δθ = -ε·tan(θ_B),
    scaling linearly with strain,
  * instrument convolution only broadens (conserves intensity, does not move
    the peak, does not invent structure).

The frozen-notebook regression (bit-for-bit match to the archival
thermo-elastic-gaas calculator) is enforced in tests/test_xrd_acceptance.py
against tests/data/gaas004_golden.npz.

Writes docs/physics_acceptance.json and, unless --no-figure, three diagnostic
panels to docs/images/xrd_acceptance.png. Exits nonzero on any failure.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from xrd_strain.acceptance import (
    kinematic_bragg_deg,
    perfect_crystal_curve,
    run_acceptance_suite,
    strained_layer_curve,
    tan_bragg,
)
from xrd_strain.config import DetectorConfig
from xrd_strain.detector.gaussian import apply_gaussian_instrument
from xrd_strain.detector.resolution import apply_detector_resolution

REPO = Path(__file__).resolve().parents[1]
REPORT = REPO / "docs" / "physics_acceptance.json"
FIGURE = REPO / "docs" / "images" / "xrd_acceptance.png"


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def _write_figure(path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    th_b = kinematic_bragg_deg()

    th_pc, i_pc = perfect_crystal_curve()
    th_l1, i_l1 = strained_layer_curve(1e-3)
    th_l2, i_l2 = strained_layer_curve(2e-3)
    th_ins, raw = strained_layer_curve(1.2e-3, n_cap_layers=60, n_angles=900)
    rad = np.radians(th_ins)
    aps = apply_gaussian_instrument(rad, raw, fwhm_arcsec=12.0)
    emp = apply_detector_resolution(rad, DetectorConfig().as_tuple(), raw)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

    ax = axes[0]
    ax.plot((th_pc - th_b) * 3600, i_pc, color="#1f77b4")
    ax.axvline(0, color="0.6", ls="--", lw=0.8, label="kinematic Bragg")
    ax.set_title("1. Perfect crystal (unstrained)")
    ax.set_xlabel("θ − θ_B (arcsec)")
    ax.set_ylabel("reflectivity")
    ax.legend(fontsize=8)

    ax = axes[1]
    for th, i, e, c in [
        (th_l1, i_l1, 1e-3, "#1f77b4"),
        (th_l2, i_l2, 2e-3, "#d62728"),
    ]:
        ax.semilogy((th - th_b) * 3600, i, color=c, lw=1.0, label=f"ε = {e:.0e}")
        ax.axvline(
            -np.degrees(e * tan_bragg()) * 3600, color=c, ls=":", lw=0.9
        )
    ax.set_title("2. Strained layer: peak at −ε·tan(θ_B)")
    ax.set_xlabel("θ − θ_B (arcsec)")
    ax.set_ylabel("intensity (log)")
    ax.legend(fontsize=8)

    ax = axes[2]
    ax.semilogy((th_ins - th_b) * 3600, raw, color="0.4", lw=1.0, label="none")
    ax.semilogy((th_ins - th_b) * 3600, aps, color="#1f77b4", lw=1.0, label="aps_7idc 12″")
    ax.semilogy((th_ins - th_b) * 3600, emp, color="#d62728", lw=1.0, label="empirical")
    ax.set_title("3. Instrument convolution (broaden only)")
    ax.set_xlabel("θ − θ_B (arcsec)")
    ax.set_ylabel("intensity (log)")
    ax.legend(fontsize=8)

    fig.suptitle("GaAs (004) XRD physics acceptance", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--no-figure", action="store_true", help="skip writing the diagnostic figure"
    )
    args = parser.parse_args()

    report = run_acceptance_suite()

    print("GaAs (004) XRD physics acceptance suite\n")
    for check in report.checks:
        status = "PASS" if check.passed else "FAIL"
        print(f"  [{status}] {check.name}: {check.detail}")

    payload = _jsonable(report.to_dict())
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"\n{payload['n_passed']}/{payload['n_total']} checks passed. "
        f"Report saved to {REPORT}"
    )

    if not args.no_figure:
        _write_figure(FIGURE)
        print(f"Figure saved to {FIGURE}")

    print("Acceptance", "PASSED" if report.passed else "FAILED")
    return 0 if report.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
