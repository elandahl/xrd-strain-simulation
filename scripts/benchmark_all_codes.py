#!/usr/bin/env python3
"""Combined three-code benchmark figure for all four APS substrates.

One panel per substrate (GaAs, Si, Ge, InSb, 004, 10 keV, sigma polarization,
symmetric Bragg).  Each panel overlays three independent dynamical-diffraction
codes for three synthetic strain profiles:

* this package's production calculator (solid, drawn first / underneath);
* xrayutilities DynamicalModel (dashed, on top);
* Stepanov GID_sl/X0h (open circle markers, on top).

Colour encodes the strain profile; line style / marker encodes the code.  No
angular or temporal instrument response is applied -- this compares diffraction
physics only.  GID_sl reference curves live in
``tests/data/gid_sl_combined/<material>_<case>.dat`` (all on one -250..+50
arcsec grid); ours and xrayutilities are evaluated on the same scan about each
code's own kinematic Bragg angle.

Outputs ``docs/images/benchmark_all_codes.png``.
"""

from __future__ import annotations

import copy
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xrayutilities as xu
from matplotlib.lines import Line2D
from xrayutilities.simpack import DynamicalModel, Layer, LayerStack

from xrd_strain.crystals.base import get_crystal
from xrd_strain.crystals.gaas_004_10kev_300k import gaas_004_10kev_300k_constants
from xrd_strain.crystals.ge_004_10kev_300k import ge_004_10kev_300k_constants
from xrd_strain.crystals.insb_004_10kev_300k import (
    insb_004_10kev_300k_constants,
)
from xrd_strain.crystals.si_004_10kev_300k import si_004_10kev_300k_constants

REPO = Path(__file__).resolve().parents[1]
GID_DIR = REPO / "tests" / "data" / "gid_sl_combined"
FIGURE = REPO / "docs" / "images" / "benchmark_all_codes.png"

E_EV = 10000.0
HKL = (0, 0, 4)
DZ_ANGSTROM = 26.7
N_CELLS = 3000
EPS = 1e-6

MATERIALS = [
    ("GaAs", "gaas", xu.materials.GaAs, gaas_004_10kev_300k_constants(), "gaas_004_10kev"),
    ("Si", "si", xu.materials.Si, si_004_10kev_300k_constants(), "si_004_10kev"),
    ("Ge", "ge", xu.materials.Ge, ge_004_10kev_300k_constants(), "ge_004_10kev"),
    ("InSb", "insb", xu.materials.InSb, insb_004_10kev_300k_constants(), "insb_004_10kev"),
]

# (label, layers as (n_cells, strain) from the surface down, colour)
CASES = [
    ("perfect crystal", [], "tab:blue"),
    ("uniform 267 nm, +1e-3", [(100, 1e-3)], "tab:orange"),
    ("two-step +2e-3 / +1e-3", [(40, 2e-3), (60, 1e-3)], "tab:green"),
]


def kinematic_bragg_deg(a: float, wavelength: float) -> float:
    return float(np.degrees(np.arcsin(wavelength / (2.0 * (a / 4.0)))))


def build_strain(layers: list[tuple[int, float]]) -> np.ndarray:
    strain = np.zeros(N_CELLS)
    start = 0
    for n_cells, value in layers:
        strain[start : start + n_cells] = value
        start += n_cells
    if strain[0] == 0.0:
        strain[0] = EPS
    return strain


def xu_curve(material, layers, scan) -> np.ndarray:
    theta_b = kinematic_bragg_deg(material.a, xu.utilities.en2lam(E_EV))
    ai = theta_b + scan / 3600.0
    stack = [Layer(material, 1e7)]
    for n_cells, eps in reversed(layers):
        mat = copy.deepcopy(material)
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
    return model.simulate(ai, hkl=HKL)


def our_curve(calculator, constants, layers, scan) -> np.ndarray:
    theta_b = kinematic_bragg_deg(
        constants.lattice_angstrom, constants.wavelength_angstrom
    )
    th = theta_b + scan / 3600.0
    return get_crystal(calculator).compute_intensity(
        th, build_strain(layers), DZ_ANGSTROM, EPS
    )


def main() -> None:
    fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharex=True)
    axes = axes.ravel()

    for ax, (name, key, material, constants, calculator) in zip(axes, MATERIALS):
        for label, layers, colour in CASES:
            scan, gid = np.loadtxt(GID_DIR / f"{key}_{label_key(label)}.dat").T
            ours = our_curve(calculator, constants, layers, scan)
            xu_int = xu_curve(material, layers, scan)

            # ours: solid, drawn first (underneath)
            ax.semilogy(scan, ours, color=colour, lw=2.2, alpha=0.9, zorder=1)
            # xrayutilities: dashed, on top
            ax.semilogy(
                scan, xu_int, color="k", lw=1.0, ls="--", alpha=0.85, zorder=2
            )
            # GID_sl: open circles, on top and clearly visible
            ax.semilogy(
                scan[::45],
                gid[::45],
                color=colour,
                marker="o",
                mfc="none",
                mec="k",
                ms=5,
                lw=0,
                zorder=3,
            )
        ax.set_title(f"{name} (004)")
        ax.set_ylim(1e-6, 2)
        ax.set_xlim(-250, 50)
        ax.grid(alpha=0.3)
        ax.set_ylabel("reflectivity")
        if name in ("Ge", "InSb"):
            ax.set_xlabel(r"$\theta - \theta_B^{\mathrm{kin}}$ (arcsec)")

    colour_handles = [
        Line2D([0], [0], color=colour, lw=2.4, label=label)
        for label, _, colour in CASES
    ]
    code_handles = [
        Line2D([0], [0], color="0.3", lw=2.2, ls="-", label="this package (solid)"),
        Line2D([0], [0], color="k", lw=1.0, ls="--", label="xrayutilities (dashed)"),
        Line2D(
            [0], [0], color="0.3", marker="o", mfc="none", mec="k", ms=6, lw=0,
            label="GID_sl / X0h (circles)",
        ),
    ]
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.9))
    leg1 = fig.legend(
        handles=colour_handles,
        title="strain profile",
        loc="upper left",
        bbox_to_anchor=(0.09, 0.955),
        ncol=3,
        frameon=True,
    )
    fig.add_artist(leg1)
    fig.legend(
        handles=code_handles,
        title="code",
        loc="upper right",
        bbox_to_anchor=(0.92, 0.955),
        ncol=3,
        frameon=True,
    )
    fig.suptitle(
        "Three-code dynamical-diffraction benchmark (004, 10 keV, "
        "sigma, no instrument response)",
        fontsize=13,
        y=0.995,
    )
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE, dpi=150)
    plt.close(fig)
    print(f"figure: {FIGURE}")


def label_key(label: str) -> str:
    return {
        "perfect crystal": "perfect",
        "uniform 267 nm, +1e-3": "uniform_layer",
        "two-step +2e-3 / +1e-3": "two_step",
    }[label]


if __name__ == "__main__":
    main()
