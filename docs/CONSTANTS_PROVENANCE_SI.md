# Si (004) @ 10 keV — constants provenance

Production calculator: `si_004_10kev` (non-default; GaAs remains the package
default). Geometry matches Sci. Rep. Fig. 2: Si (004), 10 keV, σ, symmetric
Bragg, 300 K.

| Quantity | Code value | Source |
|----------|------------|--------|
| \(a(300\,\mathrm{K})\) | 5.4309 Å | Standard RT Si lattice |
| \(\lambda\) | 1.239841984 Å | Same CODATA \(hc/E\) as GaAs |
| \(r_e\) | 2.8179403205×10⁻⁵ Å | CODATA |
| \(Z\) | 14 | — |
| \(f_0(Q_{004})\) | 7.2988 | Waasmaier & Kirfel (1995) 5-Gaussian |
| \(f'\) | +0.1940 | Henke/CXRO @ 10 keV (\(f_1-Z\)) |
| \(f''\) | 0.2136 | Henke/CXRO @ 10 keV (\(f_2\)) |
| \(B\) | 0.4613 Å² | Sears & Shelley, Acta Cryst. A47, 441 (1991) |

Structure factors (diamond, \(h+k+l=4n\)):

\[
F_0 = 8(Z+f'+if''),\qquad
F_h = 8\bigl(f_0(Q_h)+f'+if''\bigr)\,e^{-Bs^2}.
\]

## External check (GID_sl, queried 2026-07-18)

Perfect crystal (`code=Silicon`, X0h database, σ):

| Metric | GID_sl | Production |
|--------|--------|------------|
| Peak reflectivity | 0.97596 | 0.9794 |
| Darwin FWHM | 2.670″ | 2.604″ |
| Peak position | +1.325″ | +1.325″ |

Uniform strained layer (267 nm, \(da/a=+10^{-3}\)): log₁₀ RMS ≈ 0.021,
correlation ≈ 0.99999.

See `tests/test_si_004_acceptance.py` and cached `tests/data/gid_sl_si_*`.
