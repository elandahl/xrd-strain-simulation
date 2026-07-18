"""Physics acceptance diagnostics for the GaAs (004) dynamical XRD calculator.

These are *internal, physics-based* checks (no experimental data) that should be
green before external benchmarking. They mirror the strain repo's
``strain_wave.acceptance`` suite and back up the claims made in the README and
docs/INSTRUMENTS.md with quantitative, reproducible numbers:

1. Unstrained perfect-crystal rocking curve
   - the dynamical peak sits at the kinematic Bragg angle (small refraction
     shift), peak reflectivity is O(1), and the Darwin FWHM is a few arcsec.

2. Known analytic strain -> rocking curve
   - a uniformly strained surface layer produces a second peak displaced from
     the substrate by the kinematic law Δθ = -ε·tan(θ_B), and the displacement
     scales linearly with the strain ε.

3. Instrument-convolution sanity
   - a symmetric physical instrument (``aps_7idc``) only broadens: it conserves
     integrated intensity, does not move the peak, and does not invent new
     structure; the ``empirical`` effective kernel broadens without adding
     maxima either.

4. Frozen-notebook regression
   - the copied ``xrd_slab_gaas`` reproduces a checked-in golden curve computed
     from the archival thermo-elastic-gaas code (tag paper-v1.0) to machine
     precision. (Lives in ``tests/test_xrd_acceptance.py`` against
     ``tests/data/gaas004_golden.npz``.)

``scripts/validate_xrd_physics.py`` and ``tests/test_xrd_acceptance.py`` both
build on the functions here.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np

from xrd_strain.config import DetectorConfig
from xrd_strain.crystals.gaas_004_dynamical import (
    A_GAAS,
    LAMB_0,
    xrd_slab_gaas_lowmem,
)
from xrd_strain.detector.gaussian import apply_gaussian_instrument
from xrd_strain.detector.resolution import apply_detector_resolution

# GaAs (004): out-of-plane d-spacing is a/4.
D_004 = A_GAAS / 4.0
# Default layer thickness for the checks (paper grid ~2.67 nm = 26.7 A).
DEFAULT_DZ_ANGSTROM = 26.7
DEFAULT_EPS = 1e-6


def kinematic_bragg_deg() -> float:
    """Kinematic Bragg angle (deg) for GaAs (004) at the calculator's energy."""
    return float(np.degrees(np.arcsin(LAMB_0 / (2.0 * D_004))))


def tan_bragg() -> float:
    return float(np.tan(np.radians(kinematic_bragg_deg())))


def fwhm_arcsec(th_deg: np.ndarray, intensity: np.ndarray) -> float:
    """Full width at half maximum (arcsec) via linear half-max interpolation."""
    y = np.asarray(intensity, dtype=float)
    half = y.max() / 2.0
    above = np.where(y >= half)[0]
    if above.size < 2:
        return float("nan")
    lo, hi = above[0], above[-1]

    def _cross(i_out, i_in):
        if i_out == i_in:
            return th_deg[i_in]
        frac = (half - y[i_out]) / (y[i_in] - y[i_out])
        return th_deg[i_out] + frac * (th_deg[i_in] - th_deg[i_out])

    left = _cross(lo - 1, lo) if lo > 0 else th_deg[lo]
    right = _cross(hi + 1, hi) if hi < len(th_deg) - 1 else th_deg[hi]
    return float((right - left) * 3600.0)


def count_maxima(intensity: np.ndarray, rel_height: float = 0.05) -> int:
    """Number of interior local maxima above ``rel_height`` * global max."""
    y = np.asarray(intensity, dtype=float)
    thresh = rel_height * y.max()
    interior = y[1:-1]
    is_max = (interior > y[:-2]) & (interior > y[2:]) & (interior > thresh)
    return int(np.sum(is_max))


def centroid_deg(th_deg: np.ndarray, intensity: np.ndarray, rel: float = 0.5) -> float:
    """Intensity-weighted centroid (deg) over points above ``rel`` * max."""
    y = np.asarray(intensity, dtype=float)
    mask = y >= y.max() * rel
    return float(np.sum(th_deg[mask] * y[mask]) / np.sum(y[mask]))


def perfect_crystal_curve(
    n_layers: int = 1500,
    n_angles: int = 600,
    half_window_deg: float = 0.03,
    dz_angstrom: float = DEFAULT_DZ_ANGSTROM,
    eps: float = DEFAULT_EPS,
):
    """Rocking curve of an essentially unstrained thick GaAs slab.

    A single negligible top layer (ε = 2·eps) triggers the calculator's
    strained-cap bookkeeping so the thick zero-strain block below is evaluated
    as one substrate propagation — i.e. the perfect-crystal response.
    """
    th = kinematic_bragg_deg()
    th_deg = np.linspace(th - half_window_deg, th + half_window_deg, n_angles)
    strain = np.zeros(n_layers)
    strain[0] = 2.0 * eps
    intensity = xrd_slab_gaas_lowmem(th_deg, strain, dz_angstrom, eps)
    return th_deg, intensity


def strained_layer_curve(
    strain_eps: float,
    n_cap_layers: int = 100,
    n_layers: int = 1500,
    n_angles: int = 700,
    dz_angstrom: float = DEFAULT_DZ_ANGSTROM,
    eps: float = DEFAULT_EPS,
):
    """Rocking curve of a uniformly strained cap on an unstrained substrate."""
    th_b = kinematic_bragg_deg()
    predicted_shift_deg = -np.degrees(strain_eps * tan_bragg())
    th_deg = np.linspace(
        th_b + predicted_shift_deg - 0.02, th_b + 0.01, n_angles
    )
    strain = np.zeros(n_layers)
    strain[0:n_cap_layers] = strain_eps
    intensity = xrd_slab_gaas_lowmem(th_deg, strain, dz_angstrom, eps)
    return th_deg, intensity


def measure_layer_shift_arcsec(th_deg, intensity, guard_arcsec: float = 8.0) -> float:
    """Signed shift (arcsec) of the strained-layer peak from the substrate peak.

    The substrate peak is the maximum within ``guard_arcsec`` of the Bragg
    angle; the layer peak is the maximum at angles at least ``guard_arcsec``
    below it (tensile out-of-plane strain shifts the layer to lower angle).
    """
    th_b = kinematic_bragg_deg()
    guard_deg = guard_arcsec / 3600.0
    sub_region = np.abs(th_deg - th_b) < guard_deg
    i_sub = np.where(sub_region)[0][np.argmax(intensity[sub_region])]
    lay_region = th_deg < th_b - guard_deg
    if not np.any(lay_region):
        return float("nan")
    i_lay = np.where(lay_region)[0][np.argmax(intensity[lay_region])]
    return float((th_deg[i_lay] - th_deg[i_sub]) * 3600.0)


# --- Shared check/report containers (mirror strain_wave.acceptance) ----------
@dataclass
class Check:
    name: str
    passed: bool
    detail: str
    values: dict = field(default_factory=dict)


@dataclass
class AcceptanceReport:
    checks: list[Check]

    @property
    def passed(self) -> bool:
        return all(c.passed for c in self.checks)

    def to_dict(self) -> dict:
        return {
            "passed": self.passed,
            "n_passed": sum(c.passed for c in self.checks),
            "n_total": len(self.checks),
            "checks": [asdict(c) for c in self.checks],
        }


def check_perfect_crystal(n_layers: int = 1500, n_angles: int = 600) -> list[Check]:
    th_deg, intensity = perfect_crystal_curve(n_layers=n_layers, n_angles=n_angles)
    th_b = kinematic_bragg_deg()
    peak_deg = float(th_deg[int(np.argmax(intensity))])
    peak_offset_arcsec = (peak_deg - th_b) * 3600.0
    peak_reflectivity = float(intensity.max())
    width = fwhm_arcsec(th_deg, intensity)

    return [
        Check(
            name="perfect_crystal_peak_position",
            passed=abs(peak_offset_arcsec) < 3.0,
            detail=(
                f"peak at {peak_deg:.4f}° = kinematic Bragg {th_b:.4f}° "
                f"{peak_offset_arcsec:+.2f} arcsec (dynamical refraction shift)"
            ),
            values={
                "peak_deg": peak_deg,
                "bragg_deg": th_b,
                "offset_arcsec": peak_offset_arcsec,
            },
        ),
        Check(
            name="perfect_crystal_reflectivity",
            passed=0.5 < peak_reflectivity <= 1.0001,
            detail=f"peak reflectivity {peak_reflectivity:.3f} (thick perfect crystal)",
            values={"peak_reflectivity": peak_reflectivity},
        ),
        Check(
            name="perfect_crystal_darwin_width",
            passed=2.0 < width < 15.0,
            detail=f"Darwin FWHM {width:.2f} arcsec (GaAs 004, ~10 keV)",
            values={"fwhm_arcsec": width},
        ),
    ]


def check_strained_layer_shift() -> list[Check]:
    tan_b = tan_bragg()

    th1, i1 = strained_layer_curve(1e-3)
    shift1 = measure_layer_shift_arcsec(th1, i1)
    pred1 = -np.degrees(1e-3 * tan_b) * 3600.0

    th2, i2 = strained_layer_curve(2e-3)
    shift2 = measure_layer_shift_arcsec(th2, i2)
    pred2 = -np.degrees(2e-3 * tan_b) * 3600.0

    rel2 = abs(shift2 - pred2) / abs(pred2)
    scale = shift2 / shift1 if shift1 else float("nan")

    checks = [
        Check(
            name="strained_layer_bragg_shift",
            passed=(shift2 < 0) and rel2 < 0.10,
            detail=(
                f"ε=2e-3 layer peak shift {shift2:.1f} arcsec vs kinematic "
                f"-ε·tan(θ_B) = {pred2:.1f} arcsec ({rel2 * 100:.1f}% off)"
            ),
            values={"measured_arcsec": shift2, "predicted_arcsec": pred2, "rel": rel2},
        ),
        Check(
            name="strained_layer_shift_scales_with_strain",
            passed=abs(scale - 2.0) < 0.2,
            detail=(
                f"shift(2e-3)/shift(1e-3) = {scale:.2f} (expect ~2; "
                f"{shift1:.1f} -> {shift2:.1f} arcsec)"
            ),
            values={"shift_1e-3": shift1, "shift_2e-3": shift2, "scale": scale},
        ),
    ]
    return checks


def check_instrument_convolution() -> list[Check]:
    """Symmetric instrument only broadens; neither invents structure."""
    # A strained cap gives a substrate + layer peak (structure to preserve).
    th_deg, raw = strained_layer_curve(1.2e-3, n_cap_layers=60, n_angles=900)
    rad = np.radians(th_deg)

    aps = apply_gaussian_instrument(rad, raw, fwhm_arcsec=12.0)
    empirical = apply_detector_resolution(rad, DetectorConfig().as_tuple(), raw)

    fw_raw = fwhm_arcsec(th_deg, raw)
    fw_aps = fwhm_arcsec(th_deg, aps)
    fw_emp = fwhm_arcsec(th_deg, empirical)

    n_raw = count_maxima(raw)
    n_aps = count_maxima(aps)
    n_emp = count_maxima(empirical)

    area_raw = float(np.trapezoid(raw, rad))
    area_aps = float(np.trapezoid(aps, rad))
    area_ratio = area_aps / area_raw if area_raw else float("nan")

    centroid_shift = (centroid_deg(th_deg, aps) - centroid_deg(th_deg, raw)) * 3600.0

    return [
        Check(
            name="instrument_broadens",
            passed=fw_aps > fw_raw and fw_emp > fw_raw,
            detail=(
                f"FWHM none {fw_raw:.2f} -> aps(12″) {fw_aps:.2f} -> "
                f"empirical {fw_emp:.2f} arcsec"
            ),
            values={"fwhm_none": fw_raw, "fwhm_aps": fw_aps, "fwhm_empirical": fw_emp},
        ),
        Check(
            name="instrument_no_new_structure",
            passed=n_aps <= n_raw and n_emp <= n_raw,
            detail=(
                f"local maxima none {n_raw}, aps {n_aps}, empirical {n_emp} "
                "(smoothing must not add peaks)"
            ),
            values={"n_none": n_raw, "n_aps": n_aps, "n_empirical": n_emp},
        ),
        Check(
            name="aps_conserves_intensity",
            passed=abs(area_ratio - 1.0) < 0.02,
            detail=f"aps integrated-intensity ratio {area_ratio:.4f} (normalized kernel)",
            values={"area_ratio": area_ratio},
        ),
        Check(
            name="aps_does_not_move_peak",
            passed=abs(centroid_shift) < 1.0,
            detail=f"aps centroid shift {centroid_shift:+.2f} arcsec (symmetric kernel)",
            values={"centroid_shift_arcsec": centroid_shift},
        ),
    ]


def run_acceptance_suite() -> AcceptanceReport:
    checks: list[Check] = []
    checks.extend(check_perfect_crystal())
    checks.extend(check_strained_layer_shift())
    checks.extend(check_instrument_convolution())
    return AcceptanceReport(checks=checks)
