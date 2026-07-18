"""Simple Gaussian instrument-resolution model.

Used for the ``aps_7idc`` instrument: the Sci. Rep. 2022 experiment quotes a
0.5 millidegree (~1.8 arcsec) angular resolution at APS beamline 7ID-C.
"""

from __future__ import annotations

import numpy as np

ARCSEC_TO_RAD = np.pi / 180.0 / 3600.0
_FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def apply_gaussian_instrument(
    x_rad: np.ndarray, intensity: np.ndarray, fwhm_arcsec: float
) -> np.ndarray:
    """Convolve a rocking curve with a normalized Gaussian of given FWHM.

    Operates on *linear* intensity. Edge effects are handled by edge-value
    padding so the curve is not artificially darkened at the scan limits.

    Note: if the kernel FWHM is comparable to or smaller than the angular
    sampling step, the convolution is close to an identity — sample the scan
    finely enough (step < sigma) for the broadening to be meaningful.
    """
    if fwhm_arcsec <= 0:
        return intensity

    sigma = fwhm_arcsec * ARCSEC_TO_RAD * _FWHM_TO_SIGMA
    dx = float(x_rad[1] - x_rad[0])

    half_width = max(1, int(np.ceil(4.0 * sigma / dx)))
    xk = np.arange(-half_width, half_width + 1) * dx
    kernel = np.exp(-(xk**2) / (2.0 * sigma**2))
    kernel /= kernel.sum()

    padded = np.concatenate(
        (
            np.full(half_width, intensity[0]),
            intensity,
            np.full(half_width, intensity[-1]),
        )
    )
    return np.convolve(padded, kernel, mode="valid")
