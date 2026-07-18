#!/usr/bin/env python3
"""Tier-3 sensitivity of GaAs (004) rocking curves to audited constants.

Runs the identical fourth-order dynamical solver with one constant group
changed at a time from the archival notebook values to modern references:

* a(300 K): Vurgaftman et al. (2001)
* f0(Q): Waasmaier & Kirfel (1995)
* f', f''(10 keV): Henke/CXRO

The external perfect-crystal target is Stepanov's X0h server for GaAs (004) at
10 keV: sigma-polarized Darwin FWHM 5.42--5.52 arcsec across its tabulated
dispersion databases (queried 2026-07-18).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from xrd_strain.acceptance import fwhm_arcsec
from xrd_strain.crystals.gaas_004_dynamical import (
    A_GAAS,
    F_AS,
    F_GA,
    F_GAAS,
    LAMB_0,
    RE,
    xrd_slab_gaas_with_constants,
)

REPO = Path(__file__).resolve().parents[1]
REPORT = REPO / "docs" / "constants_sensitivity.json"
FIGURE = REPO / "docs" / "images" / "constants_sensitivity.png"

A_LIT = 5.65325
RE_LIT = 2.8179403205e-5
LAMBDA_LIT = 1.239841984

# Values audited in scripts/audit_constants.py.
F0_CODE = {"Ga": 19.3993, "As": 20.6498}
FP_CODE = {"Ga": -2.9830, "As": -1.4384}
FPP_CODE = {"Ga": 0.670502, "As": 0.833881}
F0_LIT = {"Ga": 19.7161, "As": 21.0693}
FP_LIT = {"Ga": -2.9739, "As": -1.6047}
FPP_LIT = {"Ga": 0.5364, "As": 0.7101}
# Room-temperature harmonic Debye-Waller B values (Angstrom^2):
# Stevenson, Acta Cryst. A50, 621-632 (1994).
B_300K = {"Ga": 0.622, "As": 0.483}

# Independent X0h result, all supported dispersion databases.
X0H_SIGMA_FWHM_RANGE_ARCSEC = (5.4238, 5.5212)
# X0h/GID_sl semi-infinite perfect-crystal curve using its recommended X0h
# dispersion database (same 10 keV, 004, sigma query).
X0H_SIGMA_PEAK_REFLECTIVITY = 0.97282552


@dataclass(frozen=True)
class ConstantsCase:
    name: str
    a: float = A_GAAS
    re: float = RE
    wavelength: float = LAMB_0
    f0_ga: float = F0_CODE["Ga"]
    f0_as: float = F0_CODE["As"]
    fp_ga: float = FP_CODE["Ga"]
    fp_as: float = FP_CODE["As"]
    fpp_ga: float = FPP_CODE["Ga"]
    fpp_as: float = FPP_CODE["As"]
    b_ga: float = 0.0
    b_as: float = 0.0

    @property
    def structure_factors(self) -> np.ndarray:
        f_ga = self.f0_ga + self.fp_ga + 1j * self.fpp_ga
        f_as = self.f0_as + self.fp_as + 1j * self.fpp_as
        s = 1.0 / (2.0 * (self.a / 4.0))
        dw_ga = np.exp(-self.b_ga * s**2)
        dw_as = np.exp(-self.b_as * s**2)
        f_004 = 4.0 * (f_ga * dw_ga + f_as * dw_as)
        return np.array([f_004, f_004], dtype=np.complex128)

    @property
    def bragg_deg(self) -> float:
        return float(
            np.degrees(np.arcsin(self.wavelength / (2.0 * (self.a / 4.0))))
        )


CASES = [
    ConstantsCase("archival"),
    ConstantsCase("a(300 K) only", a=A_LIT),
    ConstantsCase("f0 only", f0_ga=F0_LIT["Ga"], f0_as=F0_LIT["As"]),
    ConstantsCase("f' only", fp_ga=FP_LIT["Ga"], fp_as=FP_LIT["As"]),
    ConstantsCase("f'' only", fpp_ga=FPP_LIT["Ga"], fpp_as=FPP_LIT["As"]),
    ConstantsCase("Debye-Waller only", b_ga=B_300K["Ga"], b_as=B_300K["As"]),
    ConstantsCase("re + hc only", re=RE_LIT, wavelength=LAMBDA_LIT),
    ConstantsCase(
        "all audited literature",
        a=A_LIT,
        re=RE_LIT,
        wavelength=LAMBDA_LIT,
        f0_ga=F0_LIT["Ga"],
        f0_as=F0_LIT["As"],
        fp_ga=FP_LIT["Ga"],
        fp_as=FP_LIT["As"],
        fpp_ga=FPP_LIT["Ga"],
        fpp_as=FPP_LIT["As"],
    ),
    ConstantsCase(
        "all literature + DW",
        a=A_LIT,
        re=RE_LIT,
        wavelength=LAMBDA_LIT,
        f0_ga=F0_LIT["Ga"],
        f0_as=F0_LIT["As"],
        fp_ga=FP_LIT["Ga"],
        fp_as=FP_LIT["As"],
        fpp_ga=FPP_LIT["Ga"],
        fpp_as=FPP_LIT["As"],
        b_ga=B_300K["Ga"],
        b_as=B_300K["As"],
    ),
]


def perfect_curve(case: ConstantsCase, n_angles: int = 700):
    th = np.linspace(case.bragg_deg - 0.03, case.bragg_deg + 0.03, n_angles)
    strain = np.zeros(1500)
    strain[0] = 2e-6
    intensity = xrd_slab_gaas_with_constants(
        th,
        strain,
        26.7,
        1e-6,
        a_gaas=case.a,
        f_gaas=case.structure_factors,
        re=case.re,
        wavelength=case.wavelength,
    )
    return th, intensity


def layer_shift(case: ConstantsCase, strain_value: float = 2e-3) -> float:
    predicted_deg = -np.degrees(strain_value * np.tan(np.radians(case.bragg_deg)))
    th = np.linspace(
        case.bragg_deg + predicted_deg - 0.02,
        case.bragg_deg + 0.01,
        700,
    )
    strain = np.zeros(1500)
    strain[:100] = strain_value
    intensity = xrd_slab_gaas_with_constants(
        th,
        strain,
        26.7,
        1e-6,
        a_gaas=case.a,
        f_gaas=case.structure_factors,
        re=case.re,
        wavelength=case.wavelength,
    )
    guard = 8.0 / 3600.0
    sub = np.abs(th - case.bragg_deg) < guard
    layer = th < case.bragg_deg - guard
    i_sub = np.where(sub)[0][np.argmax(intensity[sub])]
    i_layer = np.where(layer)[0][np.argmax(intensity[layer])]
    return float((th[i_layer] - th[i_sub]) * 3600.0)


def analytic_darwin_width_arcsec(case: ConstantsCase) -> float:
    """Nonabsorbing symmetric-Bragg sigma width, 2|chi_h|/sin(2 theta_B)."""
    f_abs = abs(case.structure_factors[1])
    chi_h = case.re * case.wavelength**2 * f_abs / (np.pi * case.a**3)
    theta = np.radians(case.bragg_deg)
    width_rad = 2.0 * chi_h / np.sin(2.0 * theta)
    return float(np.degrees(width_rad) * 3600.0)


def _write_figure(curves, rows):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    ax = axes[0]
    colors = plt.cm.viridis(np.linspace(0.05, 0.9, len(curves)))
    for (case, th, intensity), color in zip(curves, colors):
        x = (th - case.bragg_deg) * 3600.0
        ax.plot(x, intensity, lw=1.1, color=color, label=case.name)
    ax.set_xlim(-20, 20)
    ax.set_xlabel(r"$\theta-\theta_B$ (arcsec)")
    ax.set_ylabel("reflectivity")
    ax.set_title("Perfect-crystal Darwin curves")
    ax.legend(fontsize=7)

    ax = axes[1]
    names = [row["name"] for row in rows]
    widths = [row["fwhm_arcsec"] for row in rows]
    y = np.arange(len(names))
    ax.barh(y, widths, color=colors)
    ax.axvspan(
        *X0H_SIGMA_FWHM_RANGE_ARCSEC,
        color="#d62728",
        alpha=0.18,
        label="X0h σ range",
    )
    ax.set_yticks(y, names, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Darwin FWHM (arcsec)")
    ax.set_title("One-at-a-time constants sensitivity")
    ax.legend(fontsize=8)

    fig.suptitle("GaAs (004), 10 keV — constants sensitivity")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE, dpi=160)
    plt.close(fig)


def main() -> int:
    rows = []
    curves = []
    for case in CASES:
        th, intensity = perfect_curve(case)
        curves.append((case, th, intensity))
        peak_i = int(np.argmax(intensity))
        row = {
            "name": case.name,
            "bragg_deg": case.bragg_deg,
            "peak_deg": float(th[peak_i]),
            "peak_offset_arcsec": float(
                (th[peak_i] - case.bragg_deg) * 3600.0
            ),
            "peak_reflectivity": float(intensity[peak_i]),
            "fwhm_arcsec": fwhm_arcsec(th, intensity),
            "analytic_nonabsorbing_fwhm_arcsec": analytic_darwin_width_arcsec(
                case
            ),
            "layer_shift_eps_2e3_arcsec": layer_shift(case),
            "structure_factor_abs": float(abs(case.structure_factors[1])),
        }
        rows.append(row)

    baseline = rows[0]
    for row in rows:
        row["delta_fwhm_percent_vs_archival"] = (
            100.0 * (row["fwhm_arcsec"] / baseline["fwhm_arcsec"] - 1.0)
        )
        row["delta_peak_percent_vs_archival"] = (
            100.0
            * (row["peak_reflectivity"] / baseline["peak_reflectivity"] - 1.0)
        )
        row["delta_fwhm_percent_vs_x0h_midpoint"] = 100.0 * (
            row["fwhm_arcsec"]
            / (sum(X0H_SIGMA_FWHM_RANGE_ARCSEC) / 2.0)
            - 1.0
        )
        row["delta_peak_percent_vs_x0h"] = 100.0 * (
            row["peak_reflectivity"] / X0H_SIGMA_PEAK_REFLECTIVITY - 1.0
        )

    report = {
        "case": "GaAs (004), 10 keV, sigma polarization",
        "external_reference": {
            "name": "Stepanov X0h",
            "url": "https://x-server.gmca.aps.anl.gov/x0h.html",
            "queried": "2026-07-18",
            "sigma_darwin_fwhm_arcsec_range": list(
                X0H_SIGMA_FWHM_RANGE_ARCSEC
            ),
            "sigma_peak_reflectivity_x0h_database": X0H_SIGMA_PEAK_REFLECTIVITY,
            "database_options": "X0h, Henke, Brennan-Cowan, Windt, Chantler",
        },
        "cases": rows,
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n")
    _write_figure(curves, rows)

    print("GaAs (004) @ 10 keV — Tier-3 constants sensitivity\n")
    print(
        f"{'case':26s} {'FWHM(\")':>9s} {'Δwidth':>9s} "
        f"{'Rpeak':>9s} {'Δpeak':>9s} {'layer Δθ(\")':>12s}"
    )
    for row in rows:
        print(
            f"{row['name']:26s} {row['fwhm_arcsec']:9.3f} "
            f"{row['delta_fwhm_percent_vs_archival']:+8.2f}% "
            f"{row['peak_reflectivity']:9.5f} "
            f"{row['delta_peak_percent_vs_archival']:+8.2f}% "
            f"{row['layer_shift_eps_2e3_arcsec']:12.2f}"
        )
    print(
        f"\nExternal X0h σ Darwin FWHM: "
        f"{X0H_SIGMA_FWHM_RANGE_ARCSEC[0]:.3f}–"
        f"{X0H_SIGMA_FWHM_RANGE_ARCSEC[1]:.3f} arcsec"
    )
    print(
        f"External X0h/GID_sl σ peak reflectivity: "
        f"{X0H_SIGMA_PEAK_REFLECTIVITY:.5f}"
    )
    print(f"Report saved to {REPORT}")
    print(f"Figure saved to {FIGURE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
