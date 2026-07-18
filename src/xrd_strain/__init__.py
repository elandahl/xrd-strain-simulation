"""XRD simulation from strain profiles."""

from xrd_strain import crystals as _crystals  # noqa: F401 — registers default crystals
from xrd_strain.config import DetectorConfig, XrdConfig
from xrd_strain.io import load_strain_profile
from xrd_strain.models import XrdResult
from xrd_strain.pipeline import run_xrd
from xrd_strain.temporal import (
    APS_24BUNCH_FWHM_PS,
    gaussian_delay_weights,
    run_xrd_delay_averaged,
)

__all__ = [
    "APS_24BUNCH_FWHM_PS",
    "DetectorConfig",
    "XrdConfig",
    "XrdResult",
    "gaussian_delay_weights",
    "load_strain_profile",
    "run_xrd",
    "run_xrd_delay_averaged",
]
