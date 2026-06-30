"""Dynamical X-ray diffraction simulation for strained GaAs."""

import numpy as np
from numba import njit
from numpy.linalg import inv

# Classical electron radius
RE = 2.818e-5

# GaAs lattice and scattering properties
A_GAAS = 5.65
H_GAAS = 8 * np.pi / A_GAAS * np.array([0.0, 0.0, -1.0])
F_AS = 20.6498 - 1.4384 + 1j * 0.833881
F_GA = 19.3993 - 2.9830 + 1j * 0.670502
F_GAAS = np.zeros(2, dtype=np.complex128)
F_GAAS[0] = 4 * F_AS + 4 * F_GA
F_GAAS[1] = 4 * F_AS + 4 * F_GA

SURFACE_NORMAL = np.array([0.0, 0.0, 1.0])

# Incident photon properties (10 keV)
E_0 = 10.0
LAMB_0 = 12.3984 / E_0
K_0 = 2 * np.pi / LAMB_0

M_TO_ANGSTROM = 1e10


@njit
def sort_complex_imag(compl):
    """Sort complex numbers by imaginary part in descending order."""
    for i in range(0, len(compl)):
        for j in range(i + 1, len(compl)):
            if compl[i].imag < compl[j].imag:
                temp = compl.copy()[i]
                compl[i] = compl.copy()[j]
                compl[j] = temp

    return compl.copy()


@njit
def solve4thOrder(az, daz, a_sub, F, psi, gamma_0, gamma_h):
    psi_n = psi * (1 - daz / az)

    gam = RE * (LAMB_0**2) / (np.pi * (a_sub**2) * az)

    chi = np.zeros((2), dtype=np.complex128)

    chi[0] = -gam * F[0].real + gam * np.abs(F[0].imag) * 1j
    chi[1] = -gam * F[1].real + gam * np.abs(F[1].imag) * 1j

    a = 1.0 + np.zeros(len(gamma_0), dtype=np.complex128)
    b = 2 * psi_n * np.ones(len(gamma_0), dtype=np.complex128)
    c = psi_n**2 - gamma_0**2 - gamma_h**2 - 2 * chi[0]
    d = -2 * psi_n * (gamma_0**2 + chi[0])
    e = (
        (gamma_0**2) * (gamma_h**2 - psi_n**2)
        + (gamma_0**2 + gamma_h**2 - psi_n**2) * chi[0]
        + chi[0] ** 2
        - chi[1] ** 2
    )

    u = np.zeros((len(gamma_0), 4), dtype=np.complex128)
    v = np.zeros((len(gamma_0), 4), dtype=np.complex128)
    w = np.zeros((len(gamma_0), 4), dtype=np.complex128)

    for i in range(0, len(gamma_0)):
        u[i] = sort_complex_imag(
            np.roots(np.array([a[i], b[i], c[i], d[i], e[i]], dtype=np.complex128))
        ).copy()

    for i in range(0, len(gamma_0)):
        v[i] = ((u[i] ** 2 - gamma_0[i] ** 2 - chi[0]) / chi[1]).copy()
        w[i] = (v[i] * (u[i] + psi_n)).copy()

    return u.copy(), v.copy(), w.copy()


@njit
def xrd_slab_gaas(th_0, strain, dz, eps):
    rad_0 = th_0 * np.pi / 180

    k_0_vec = np.zeros((3, len(rad_0)), dtype=np.float64)
    k_0_vec[0] = K_0 * np.cos(rad_0)
    k_0_vec[1] = np.zeros(len(rad_0), dtype=np.float64)
    k_0_vec[2] = K_0 * np.sin(rad_0)

    gamma_0 = np.sin(rad_0)
    psi = H_GAAS @ SURFACE_NORMAL / K_0
    alpha = (2 * k_0_vec.T @ H_GAAS + H_GAAS @ H_GAAS) / (K_0**2)
    gamma_h = np.sqrt((gamma_0 + psi) ** 2 - alpha)

    Ni = np.zeros(len(strain))
    Nf = np.zeros(len(strain))
    dN = np.zeros(len(strain))

    j = 0
    ki = 0
    kf = 0
    for i in range(0, int(len(strain))):
        j_pre = j

        if np.abs(strain[i]) < eps and ki == 0:
            Ni[j] = i
            ki += 1
        elif np.abs(strain[i]) < eps and ki != 0:
            ki += 1
        elif np.abs(strain[i]) > eps and ki != 0 and kf == 0:
            Nf[j] = i
            dN[j] = Nf[j] - Ni[j]
            kf += 1
            j += 1

        if j_pre != j:
            ki = 0
            kf = 0

        if len(Ni[Ni != 0]) != len(Nf[Nf != 0]):
            Nf[j] = len(strain)
            dN[j] = Nf[j] - Ni[j]

    N_layer = int(len(strain))
    N_memory = N_layer - int(dN.sum()) + len(Ni[Ni != 0]) + 2

    S = np.zeros((N_memory + 1, len(gamma_0), 4, 4), dtype=np.complex128)
    S_inv = np.zeros((N_memory + 1, len(gamma_0), 4, 4), dtype=np.complex128)
    A = np.zeros((N_memory + 1, len(gamma_0), 4, 4), dtype=np.complex128)

    W = np.zeros((len(gamma_0), 4, 4), dtype=np.complex128)
    d = np.zeros(len(gamma_0), dtype=np.complex128)
    r = np.zeros(len(gamma_0), dtype=np.complex128)

    u = np.zeros((N_memory + 1, len(gamma_0), 4), dtype=np.complex128)
    v = np.zeros((N_memory + 1, len(gamma_0), 4), dtype=np.complex128)
    w = np.zeros((N_memory + 1, len(gamma_0), 4), dtype=np.complex128)

    l = 0
    m = 1

    for i in range(0, N_layer):
        if i <= Ni[l]:
            da = strain[i] * A_GAAS
            u[m], v[m], w[m] = solve4thOrder(
                A_GAAS, da, A_GAAS, F_GAAS, psi, gamma_0, gamma_h
            )
            m += 1
        elif i > Ni[l] and i < Nf[l]:
            pass
        elif i == Nf[l]:
            da = strain[i] * A_GAAS
            u[m], v[m], w[m] = solve4thOrder(
                A_GAAS, da, A_GAAS, F_GAAS, psi, gamma_0, gamma_h
            )
            m += 1
            l += 1

    for j_idx in range(0, len(gamma_0)):
        S[0][j_idx][0] = np.array([1, 0, 1, 0])
        S[0][j_idx][1] = np.array([0, 1, 0, 1])
        S[0][j_idx][2] = np.array([gamma_0[j_idx], 0, -gamma_0[j_idx], 0])
        S[0][j_idx][3] = np.array([0, gamma_h[j_idx], 0, -gamma_h[j_idx]])
        S_inv[0][j_idx] = inv(S[0][j_idx])

        l = 0
        m = 1

        for i in range(0, N_layer):
            if i < Ni[l]:
                S[m][j_idx][0] = np.array([1, 1, 1, 1])
                S[m][j_idx][1] = v[m][j_idx].copy()
                S[m][j_idx][2] = u[m][j_idx].copy()
                S[m][j_idx][3] = w[m][j_idx].copy()
                S_inv[m][j_idx] = inv(S[m][j_idx])

                A[m][j_idx][0][0] = np.exp(-1j * u[m][j_idx][0] * K_0 * dz)
                A[m][j_idx][1][1] = np.exp(-1j * u[m][j_idx][1] * K_0 * dz)
                A[m][j_idx][2][2] = np.exp(-1j * u[m][j_idx][2] * K_0 * dz)
                A[m][j_idx][3][3] = np.exp(-1j * u[m][j_idx][3] * K_0 * dz)

                if i == 0:
                    W[j_idx] = S_inv[m - 1][j_idx] @ S[m][j_idx] @ A[m][j_idx]
                else:
                    W[j_idx] = W[j_idx] @ S_inv[m - 1][j_idx] @ S[m][j_idx] @ A[m][j_idx]

                m += 1

            elif i == Ni[l]:
                S[m][j_idx][0] = np.array([1, 1, 1, 1])
                S[m][j_idx][1] = v[m][j_idx].copy()
                S[m][j_idx][2] = u[m][j_idx].copy()
                S[m][j_idx][3] = w[m][j_idx].copy()
                S_inv[m][j_idx] = inv(S[m][j_idx])

                A[m][j_idx][0][0] = np.exp(-1j * u[m][j_idx][0] * K_0 * dz * dN[l])
                A[m][j_idx][1][1] = np.exp(-1j * u[m][j_idx][1] * K_0 * dz * dN[l])
                A[m][j_idx][2][2] = np.exp(-1j * u[m][j_idx][2] * K_0 * dz * dN[l])
                A[m][j_idx][3][3] = np.exp(-1j * u[m][j_idx][3] * K_0 * dz * dN[l])

                if i == 0:
                    W[j_idx] = S_inv[m - 1][j_idx] @ S[m][j_idx] @ A[m][j_idx]
                else:
                    W[j_idx] = W[j_idx] @ S_inv[m - 1][j_idx] @ S[m][j_idx] @ A[m][j_idx]

                m += 1

            elif i > Ni[l] and i < Nf[l]:
                pass

            elif i == Nf[l]:
                S[m][j_idx][0] = np.array([1, 1, 1, 1])
                S[m][j_idx][1] = v[m][j_idx].copy()
                S[m][j_idx][2] = u[m][j_idx].copy()
                S[m][j_idx][3] = w[m][j_idx].copy()
                S_inv[m][j_idx] = inv(S[m][j_idx])

                A[m][j_idx][0][0] = np.exp(-1j * u[m][j_idx][0] * K_0 * dz)
                A[m][j_idx][1][1] = np.exp(-1j * u[m][j_idx][1] * K_0 * dz)
                A[m][j_idx][2][2] = np.exp(-1j * u[m][j_idx][2] * K_0 * dz)
                A[m][j_idx][3][3] = np.exp(-1j * u[m][j_idx][3] * K_0 * dz)

                W[j_idx] = W[j_idx] @ S_inv[m - 1][j_idx] @ S[m][j_idx] @ A[m][j_idx]

                l += 1
                m += 1

        D1 = -W[j_idx][1][1] / (
            W[j_idx][0][1] * W[j_idx][1][0] - W[j_idx][1][1] * W[j_idx][0][0]
        )
        D2 = W[j_idx][1][0] / (
            W[j_idx][0][1] * W[j_idx][1][0] - W[j_idx][1][1] * W[j_idx][0][0]
        )
        d[j_idx] = D1 * W[j_idx][3][0] + D2 * W[j_idx][3][1]
        r[j_idx] = D1 * W[j_idx][2][0] + D2 * W[j_idx][2][1]

    intensity = np.abs(d) ** 2 * gamma_h / gamma_0
    return intensity


def _gaussian_filter(x, mean, sig):
    return np.exp(-((x - mean) ** 2) / (2 * sig**2)) / (np.sqrt(2 * np.pi) * sig)


def apply_detector_resolution(x, detector_params, y_fit):
    """Convolve simulated rocking curve with instrument resolution."""
    a1, a2, a3, cen1, cen2, cen3, cen4, p1, p2, p3, p4 = detector_params

    sig_conv1 = p1 * np.pi / 180
    sig_conv2 = p2 * np.pi / 180
    sig_conv3 = p3 * np.pi / 180
    sig_conv4 = p4 * np.pi / 180

    a4 = 1 - a1 - a2 - a3

    dx = x[1] - x[0]
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

    return y_fit[n_1 : n_1 + len(x)]


def detector_params_from_config(detector) -> tuple:
    return (
        detector.a1,
        detector.a2,
        detector.a3,
        detector.cen1,
        detector.cen2,
        detector.cen3,
        detector.cen4,
        detector.p1,
        detector.p2,
        detector.p3,
        detector.p4,
    )
