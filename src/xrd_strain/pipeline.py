"""XRD pipeline."""

import numpy as np

from xrd_strain.config import M_TO_ANGSTROM, DetectorConfig, XrdConfig
from xrd_strain.crystals.base import get_crystal
from xrd_strain.detector.gaussian import apply_gaussian_instrument
from xrd_strain.detector.resolution import apply_detector_resolution
from xrd_strain.io import StrainProfile
from xrd_strain.models import XrdResult

INSTRUMENTS = ("empirical", "aps_7idc", "none")

# Backward-compatible aliases for renamed instruments.
_INSTRUMENT_ALIASES = {"notebook": "empirical"}


def resolve_instrument(name: str) -> str:
    """Map deprecated instrument names to their current key."""
    return _INSTRUMENT_ALIASES.get(name, name)


def run_xrd(
    profile: StrainProfile,
    config: XrdConfig | None = None,
    detector: DetectorConfig | None = None,
) -> XrdResult:
    config = config or XrdConfig()
    detector = detector or DetectorConfig()

    crystal = get_crystal(config.crystal)
    if profile.substrate_material != crystal.substrate_material:
        raise ValueError(
            f"Strain profile substrate {profile.substrate_material!r} does not match "
            f"crystal calculator substrate {crystal.substrate_material!r}."
        )

    default_min, default_max = crystal.default_angle_range_deg
    angle_min = default_min if config.angle_min is None else config.angle_min
    angle_max = default_max if config.angle_max is None else config.angle_max
    th_deg = np.linspace(angle_min, angle_max, config.n_points)
    rad = th_deg * np.pi / 180

    intensity = crystal.compute_intensity(
        th_deg,
        profile.substrate_strain,
        profile.dz * M_TO_ANGSTROM,
        config.strain_eps,
    )

    instrument = resolve_instrument(config.instrument)
    if instrument == "empirical":
        intensity = apply_detector_resolution(rad, detector.as_tuple(), intensity)
    elif instrument == "aps_7idc":
        intensity = apply_gaussian_instrument(
            rad, intensity, config.instrument_fwhm_arcsec
        )
    elif instrument == "none":
        pass
    else:
        raise ValueError(
            f"Unknown instrument {config.instrument!r}. Available: {INSTRUMENTS}"
        )

    if config.log10_intensity:
        intensity = np.log10(intensity)

    return XrdResult(
        angle_deg=th_deg,
        intensity=intensity,
        crystal=config.crystal,
        substrate_material=profile.substrate_material,
        strain_model=profile.model,
        instrument=instrument,
    )
