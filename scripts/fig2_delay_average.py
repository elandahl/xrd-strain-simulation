"""Fig. 2 (Cr/Si) forward model with APS temporal resolution.

The published Sci. Rep. Fig. 2(a) rocking curve at nominal delay 0.34 ns was
measured with an x-ray bunch duration of ~90 ps FWHM (APS 24-bunch mode).
During one bunch the acoustic wavefront in Si moves v_Si * 90 ps ~ 760 nm, so
the measurement is an incoherent average of rocking curves over that delay
window.  Two consequences:

- The coherent depth (Pendelloesung-like) fringes, whose ~5 arcsec spacing is
  set by the instantaneous strained depth, scramble and cancel.
- The Brillouin sidebands survive, because their angular position is fixed by
  the acoustic pulse-train period (~200 nm), which does not change with delay.

This script:
1. runs the `paper_fig2_si` strain preset at N Gaussian-weighted delays
   spanning +/- 2 sigma about 0.34 ns (via the sibling strain-wave-simulation
   repo, using its venv python if present);
2. computes the Si (004) rocking curve for each delay with the production
   `si_004_10kev` calculator;
3. forms the weighted average, applies the aps_7idc 1.8 arcsec instrument, and
   writes `docs/images/fig2_delay_average.png` plus machine-readable metrics
   into `docs/fig2_forward_summary.json`.

Run from the repo root:  python scripts/fig2_delay_average.py
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
STRAIN_REPO = REPO.parent / "strain-wave-simulation"
DELAY_DIR = STRAIN_REPO / "results" / "paper_fig2_si" / "delays"

BUNCH_FWHM_S = 90e-12       # APS 24-bunch x-ray duration
CENTER_DELAY_S = 0.34e-9    # nominal paper delay
N_DELAYS = 9                # samples over +/- 2 sigma
ARCSEC_HALF_RANGE = 160
N_ANGLES = 641
INSTRUMENT_FWHM_ARCSEC = 1.8


def python_in(repo: Path) -> str:
    venv = repo / ".venv" / "bin" / "python"
    return str(venv) if venv.exists() else sys.executable


def run_strain_delays() -> tuple[np.ndarray, np.ndarray]:
    """Run the strain preset at each delay; return (times_s, weights)."""
    sigma = BUNCH_FWHM_S / (2 * np.sqrt(2 * np.log(2)))
    offsets = np.linspace(-2 * sigma, 2 * sigma, N_DELAYS)
    weights = np.exp(-(offsets**2) / (2 * sigma**2))
    weights /= weights.sum()
    times = CENTER_DELAY_S + offsets

    DELAY_DIR.mkdir(parents=True, exist_ok=True)
    code = (
        "import sys, numpy as np\n"
        "from pathlib import Path\n"
        "from strain_wave import get_preset, run_simulation, save_strain_profile\n"
        f"times = {list(map(float, times))}\n"
        f"out = Path({str(DELAY_DIR)!r})\n"
        "for k, t in enumerate(times):\n"
        "    config = get_preset('paper_fig2_si')\n"
        "    config.t_max = float(t)\n"
        "    r = run_simulation(config=config, verbose=False)\n"
        "    save_strain_profile(r.to_profile(), out / f'strain_{k}.npz')\n"
        "    print(f'  delay {k}: t={t*1e12:.1f} ps  reach={r.elastic_wave_reach_nm} nm')\n"
    )
    print(f"Running {N_DELAYS} strain simulations in {STRAIN_REPO.name} ...")
    subprocess.run([python_in(STRAIN_REPO), "-c", code], check=True, cwd=STRAIN_REPO)
    np.savez(
        DELAY_DIR / "meta.npz",
        t=times, w=weights, fwhm_s=BUNCH_FWHM_S, center_s=CENTER_DELAY_S,
    )
    return times, weights


def local_maxima(x, y, lo, hi, thresh=1e-4):
    m = (x > lo) & (x < hi)
    xs, ys = x[m], y[m]
    return [
        (float(xs[i]), float(ys[i]))
        for i in range(1, len(ys) - 1)
        if ys[i] > ys[i - 1] and ys[i] > ys[i + 1] and ys[i] > thresh
    ]


def main() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image

    from xrd_strain.crystals.base import get_crystal
    from xrd_strain.crystals.si_004_10kev_300k import si_004_10kev_300k_constants
    from xrd_strain.detector.gaussian import apply_gaussian_instrument
    from xrd_strain.io import load_strain_profile

    times, weights = run_strain_delays()

    c = si_004_10kev_300k_constants()
    theta_b = float(
        np.degrees(np.arcsin(c.wavelength_angstrom / (2 * (c.lattice_angstrom / 4))))
    )
    th = np.linspace(
        theta_b - ARCSEC_HALF_RANGE / 3600,
        theta_b + ARCSEC_HALF_RANGE / 3600,
        N_ANGLES,
    )
    calc = get_crystal("si_004_10kev")

    curves = []
    for k in range(len(weights)):
        profile = load_strain_profile(DELAY_DIR / f"strain_{k}.npz")
        intensity = calc.compute_intensity(
            th, profile.substrate_strain, profile.dz * 1e10, 1e-6
        )
        curves.append(intensity)
        print(f"XRD delay {k}: t={times[k]*1e12:.0f} ps done")
    curves = np.asarray(curves)
    intensity_avg = np.tensordot(weights, curves, axes=1)

    rad = th * np.pi / 180
    intensity_avg_aps = apply_gaussian_instrument(
        rad, intensity_avg, INSTRUMENT_FWHM_ARCSEC
    )
    # index of the nominal (center) delay for the single-shot reference
    k_center = int(np.argmin(np.abs(times - CENTER_DELAY_S)))
    intensity_single_aps = apply_gaussian_instrument(
        rad, curves[k_center], INSTRUMENT_FWHM_ARCSEC
    )

    x = (th - theta_b) * 3600
    out = REPO / "results" / "paper_fig2_si"
    out.mkdir(parents=True, exist_ok=True)
    np.savez(
        out / "xrd_delay_averaged.npz",
        angle_deg=th, intensity_raw=intensity_avg, intensity_aps=intensity_avg_aps,
        theta_b=theta_b, weights=weights, times=times,
    )

    pk_avg = local_maxima(x, intensity_avg_aps / intensity_avg_aps.max(), 25, 130)
    print("delay-averaged +side maxima (arcsec):", [round(p[0], 1) for p in pk_avg])

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    ax = axes[0]
    ax.semilogy(
        x, intensity_single_aps / intensity_single_aps.max(), color="0.65", lw=0.9,
        label="single delay 0.34 ns (5\u2033 depth fringes)",
    )
    ax.semilogy(
        x, intensity_avg_aps / intensity_avg_aps.max(), "r-", lw=1.5,
        label="90 ps FWHM bunch average",
    )
    for s in (-61, 61):
        ax.axvline(s, color="0.4", ls="--", lw=0.9)
    ax.text(61, 0.5, " \u00b161\u2033 (paper)", fontsize=8, color="0.3")
    ax.set_xlim(-ARCSEC_HALF_RANGE, ARCSEC_HALF_RANGE)
    ax.set_ylim(1e-3, 1.2)
    ax.set_xlabel(r"$\theta_i-\theta_B$ (arcsec)")
    ax.set_ylabel("norm. intensity")
    ax.set_title("Si (004), \u0394t = 0.34 ns, aps_7idc 1.8\u2033")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)

    ax = axes[1]
    pub = np.asarray(Image.open(REPO / "docs/images/published_fig2_scirep_2022.png"))
    ax.imshow(pub)
    ax.axis("off")
    ax.set_title("Published Sci. Rep. Fig. 2")

    fig.suptitle(
        "Fig. 2 forward model with APS 24-bunch temporal resolution "
        "(90 ps FWHM x-ray bunch)",
        fontsize=11,
    )
    fig.tight_layout()
    for p in (
        REPO / "docs/images/fig2_delay_average.png",
        STRAIN_REPO / "docs/images/fig2_delay_average.png",
        out / "fig2_delay_average.png",
    ):
        fig.savefig(p, dpi=150)
    plt.close(fig)

    summary_path = REPO / "docs" / "fig2_forward_summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    summary["temporal_resolution_model"] = {
        "bunch_fwhm_ps": BUNCH_FWHM_S * 1e12,
        "n_delay_samples": N_DELAYS,
        "delay_range_ps": [float(times[0] * 1e12), float(times[-1] * 1e12)],
        "surviving_positive_maxima_arcsec": [round(p[0], 1) for p in pk_avg],
        "note": (
            "Gaussian average over the APS 24-bunch x-ray duration removes the "
            "5-arcsec coherent depth fringes; Brillouin sidebands survive "
            "because their angle is set by the pulse-train period."
        ),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print("wrote", summary_path)


if __name__ == "__main__":
    main()
