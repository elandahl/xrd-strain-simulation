"""Temporal instrument response (delay averaging).

A time-resolved rocking curve measurement has *two* instrument responses:

- **angular**: monochromator/analyzer/detector acceptance, applied as a
  convolution of the rocking curve in theta (see ``xrd_strain.detector``);
- **temporal**: the probe duration (x-ray bunch length) plus pump-probe
  timing jitter, applied as a convolution of the *time-dependent* signal in
  delay.

The temporal response cannot be applied to a single strain snapshot: it is a
weighted incoherent average of rocking curves computed from strain profiles
at several delays around the nominal one,

    I_meas(theta; t0) = sum_k w_k * I(theta; t_k),

with w_k a normalized sampling of the temporal kernel. Since the angular
convolution is linear, the two responses commute; averaging must however
happen on *linear* intensity, before any log10.

For the APS 24-bunch mode used in the Sci. Rep. 2022 experiment the temporal
kernel is dominated by the x-ray bunch duration, ~90 ps FWHM (timing jitter
is negligible in comparison).
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from xrd_strain.config import DetectorConfig, XrdConfig
from xrd_strain.io import StrainProfile
from xrd_strain.models import XrdResult

#: X-ray bunch duration in APS 24-bunch mode (Sci. Rep. 2022 experiment).
APS_24BUNCH_FWHM_PS = 90.0

_FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def gaussian_delay_weights(
    center_s: float,
    fwhm_s: float,
    n_samples: int = 9,
    span_sigmas: float = 2.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample a Gaussian temporal kernel.

    Returns ``(times_s, weights)`` with ``n_samples`` delays spanning
    ``center_s +/- span_sigmas * sigma`` and weights normalized to sum to 1.
    With ``fwhm_s <= 0`` this degenerates to the single nominal delay.
    """
    if fwhm_s <= 0 or n_samples < 2:
        return np.array([center_s]), np.array([1.0])
    sigma = fwhm_s * _FWHM_TO_SIGMA
    offsets = np.linspace(-span_sigmas * sigma, span_sigmas * sigma, n_samples)
    weights = np.exp(-(offsets**2) / (2.0 * sigma**2))
    weights /= weights.sum()
    return center_s + offsets, weights


def run_xrd_delay_averaged(
    profiles: Sequence[StrainProfile],
    weights: np.ndarray,
    config: XrdConfig | None = None,
    detector: DetectorConfig | None = None,
) -> XrdResult:
    """Apply the temporal instrument response to a set of delay snapshots.

    ``profiles`` are strain profiles computed at the delays returned by
    :func:`gaussian_delay_weights` (in the same order as ``weights``). Each
    snapshot's rocking curve is computed with the configured *angular*
    response, then averaged on linear intensity; ``log10_intensity`` is
    applied only to the average.
    """
    from xrd_strain.pipeline import run_xrd

    config = config or XrdConfig()
    if len(profiles) != len(weights):
        raise ValueError(
            f"Got {len(profiles)} profiles but {len(weights)} weights."
        )

    snapshot_config = XrdConfig(**{**config.__dict__, "log10_intensity": False})
    result = None
    intensity = None
    for profile, w in zip(profiles, weights):
        result = run_xrd(profile, config=snapshot_config, detector=detector)
        term = w * result.intensity
        intensity = term if intensity is None else intensity + term

    if config.log10_intensity:
        intensity = np.log10(intensity)

    return XrdResult(
        angle_deg=result.angle_deg,
        intensity=intensity,
        crystal=result.crystal,
        substrate_material=result.substrate_material,
        strain_model=result.strain_model,
        instrument=result.instrument,
    )
