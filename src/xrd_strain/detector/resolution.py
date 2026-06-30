"""Detector resolution convolution."""

import numpy as np


def _gaussian_filter(x, mean, sig):
    return np.exp(-((x - mean) ** 2) / (2 * sig**2)) / (np.sqrt(2 * np.pi) * sig)


def apply_detector_resolution(x_rad: np.ndarray, detector_params: tuple, y_fit: np.ndarray) -> np.ndarray:
    a1, a2, a3, cen1, cen2, cen3, cen4, p1, p2, p3, p4 = detector_params

    sig_conv1 = p1 * np.pi / 180
    sig_conv2 = p2 * np.pi / 180
    sig_conv3 = p3 * np.pi / 180
    sig_conv4 = p4 * np.pi / 180

    a4 = 1 - a1 - a2 - a3

    dx = x_rad[1] - x_rad[0]
    n_pad = int(12 * (sig_conv1 + sig_conv2) / dx)
    n_1 = int(n_pad / 2)
    n_2 = n_1

    x_conv = np.linspace(0 - n_1 * dx, 0 + n_2 * dx, n_pad)

    c1 = cen1 * np.pi / 180
    c2 = cen2 * np.pi / 180
    c3 = cen3 * np.pi / 180
    c4 = cen4 * np.pi / 180

    y_filter = (
        a1 * _gaussian_filter(x_conv, c1, sig_conv1)
        + a2 * _gaussian_filter(x_conv, c2, sig_conv2)
        + a3 * _gaussian_filter(x_conv, c3, sig_conv3)
        + a4 * _gaussian_filter(x_conv, c4, sig_conv4)
    )

    y_fit = np.convolve(y_filter, y_fit)
    y_fit = y_fit / len(y_fit)

    return y_fit[n_1 : n_1 + len(x_rad)]
