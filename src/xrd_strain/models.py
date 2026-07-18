"""XRD result container."""

from dataclasses import dataclass

import numpy as np


@dataclass
class XrdResult:
    angle_deg: np.ndarray
    intensity: np.ndarray
    crystal: str
    substrate_material: str
    strain_model: str
    instrument: str = "aps_7idc"
