"""XRD pipeline."""

import numpy as np

from xrd_strain.config import M_TO_ANGSTROM, DetectorConfig, XrdConfig
from xrd_strain.crystals.base import get_crystal
from xrd_strain.detector.resolution import apply_detector_resolution
from xrd_strain.io import StrainProfile
from xrd_strain.models import XrdResult


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

    th_deg = np.linspace(config.angle_min, config.angle_max, config.n_points)
    rad = th_deg * np.pi / 180

    intensity = crystal.compute_intensity(
        th_deg,
        profile.substrate_strain,
        profile.dz * M_TO_ANGSTROM,
        config.strain_eps,
    )

    intensity = apply_detector_resolution(rad, detector.as_tuple(), intensity)

    if config.log10_intensity:
        intensity = np.log10(intensity)

    return XrdResult(
        angle_deg=th_deg,
        intensity=intensity,
        crystal=config.crystal,
        substrate_material=profile.substrate_material,
        strain_model=profile.model,
    )
